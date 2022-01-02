% test1 = import_png_bw("test1.png");
% [test1c_y, test1c_Cr, test1c_Cb] = import_png_ycrcb("test1c.png");

evaluate_compression_bw("test1.png") % 0x 
evaluate_compression_rgb("test1c.png")% 2x
evaluate_compression_bw("test2.png") % 1x 
evaluate_compression_rgb("test2c.png") % 2x
evaluate_compression_bw("test3.png") % 1x
evaluate_compression_rgb("test3c.png") % 0x

function evaluate_compression_bw(name)
    dat = import_png_bw(name);
    [q_prob, c] = dct_quantization(name,dat);
    fprintf("%10s:  \t%dx compression - p quant = %f\n", name, c, q_prob);
end
function evaluate_compression_rgb(name)
    [Y, Cr, Cb] = import_png_ycrcb(name);
    [q_prob, c] = dct_quantization(name,Y);
    fprintf("%10s Y:\t%dx compression - p quant = %f\n", name, c, q_prob);
    % As noted in [Chen2011], Cr and Cb are "coarsely sampled", so their 
    % "periodicity tends to be poorly characterized".
    % Basically, they don't help detect the compression.
%     [q_prob, c] = dct_quantization(name,Cr);
%     fprintf("%s Cr: %d x compression - p of quant = %f\n", name, c, q_prob);
%     [q_prob, c] = dct_quantization(name,Cb);
%     fprintf("%s Cb: %d x compression - p of quant = %f\n", name, c, q_prob);


%     q_prob = dct_quantization(Cr);
%     fprintf("Cr: p of quant = %f\n", q_prob);
%     q_prob = dct_quantization(Cb);
%     fprintf("Cb: p of quant = %f\n", q_prob);
end
function [quant_prob,c]=dct_quantization(name,dat)    
    % First step - detect any compression at all
    % Use Neelamani et al. 2006 (https://doi.org/10/bnm5nv) method for
    % color, which refines Fan et al. 2003 (https://doi.org/10/bjrkbd)
    % 1. Split into 8x8 blocks
    bs = blocks_8x8(dat);
    % 2. (2003 method only??) Filter blocks
%     filtered_bs = cell(1, length(bs));
%     filtered_i = 1;
%     for i = 1:length(bs)
%         b = bs{i};
%         % Removing potentially-truncated (those containing 0 or 255)
%         % or uniform (all elements are the same) blocks
%         if ismember(255, b) || ismember(0, b) || length(unique(b)) == 1
%             % This block is potentially-truncated, or uniform
%             % => Don't add to filtered_bs
%         else
%             filtered_bs{filtered_i} = b;
%             filtered_i = filtered_i + 1;
%         end
%     end
%     % Trim elements off
%     filtered_bs = filtered_bs(1:(filtered_i-1));
%     fprintf("Filtered %d initial blocks down to %d\n", length(bs), length(filtered_bs));
    % 3. Compute the set Di of observed DCT coefficients for each freq i
    Di = cell(8, 8); % Di = 8x8 cell array, each cell = 1xlength(filtered_bs) vector
    for j = 1:8
        for k = 1:8
            Di{j,k} = zeros(1, length(bs));
        end
    end
    for i = 1:length(bs)
        b = bs{i};
        d = round(dct2(b));
        % For unsigned 8-bit input, DCT coefficients are signed 11-bit ints
        
        % TODO - do 2003 filtering here?
        if ismember(255, b) || ismember(0, b) || length(unique(b)) == 1
            % This block is potentially-truncated, or uniform
            % => DCT could throw off results
        else
            for j = 1:8
                for k = 1:8
                    Di{j,k}(i) = d(j,k);
                end
            end
        end
    end
    % Estimate parameter lambda from observations
    lambdas = zeros(8, 8);
    probably_quantized = zeros(8,8);
%     figure;
    for j = 1:8
        for k = 1:8
            
            
            % Experiment - try doubly compressed detection?
            if j == 3 && k == 2
%             subplot(8,8,(j-1)*8+k)
%             plot(abs(fft(histcounts(Di{j,k},'BinMethod','Integer'))))
                c=evaluate_compression_level(Di{j,k});
            end
            
            Di{j,k} = abs(Di{j,k});

            
%             Di{j,k} = unique(Di{j,k});
            
            
            
            % HACK TIME
            % If this has been quantized to some level q
            % as per [Neelamani2006] the histogram of |Di| will have spikes
            % around n*q
            % => histogram of uniq(|Di|) will be strings of ones and zeros,
            % where the first string of ones = the values around Di = 0
            % and the second string of ones = the values around Di = q
            % IF we're quantized, these strings should be around 3-4
            % elements long
            uniq_di_hist = zeros(1, max(Di{j,k}));%ismember(Di{j,k}, 0:max(Di{j,k}))
            uniq_di_hist(nonzeros(Di{j,k})) = 1;
            uniq_di_hist;
            
            ls = lengths_of_one_runs(uniq_di_hist);
            avg_run_length = mean(ls);
            % BIG hack - assume all runs will be short in a quantized
            % scenario
            % Also if the first run is long it's not possible
            probably_quantized(j,k) = ls(1) < 6 && avg_run_length < 6;
            
            lambdas(j,k) = length(Di{j,k})/sum(abs(Di{j,k}));
            

            % Estimate quantization level q where q = 1-100
%             for q = 1:100
%                 
%             end
        end
    end
    sgtitle(name);
%     probably_quantized;
%     histogram(Di{4,4},'BinMethod','Integer')
    
    % JPEG compression with Q = 100 != lossless
    % https://stackoverflow.com/a/48578831/4248422
    
    quant_prob = mean(probably_quantized, 'all');
end
function c=evaluate_compression_level(di_3_2_values)
    % make histogram from abs(di_3_2_values)
    h = histcounts(abs(di_3_2_values),'BinMethod','Integer');
    % take fft of it - this will be two-sided, because h is real-valued,
    % so take only one side
    % TODO - if we do this, then surely index 0 should always be a turning
    % point? would be better to just not do it?
    y_twosided = abs(fft(h));
    y = y_twosided(ceil(length(h)/2):end);
    % the graphs have three distinct shapes for 0x, 1x, 2x compression
    % we can differentiate between them using the turning points
    ps = get_turning_points(y);
    
    figure;
    plot(y);
    hold on;
    if length(ps) > 0
        plot(ps(:, 1),ps(:, 2),'r*');
    end
    hold off;
    
    % First, if there aren't any oscillations, this is not compressed
    if length(ps) < 3
        c = 0;
    else
        % OK, we're compressed, but *how* compressed are we?
        % As observed in [Chen2011], the size of the peak magnitudes will
        % be different if compressed 2x
        % => take the differences between adjacent turning points 
        %      e.g. the "size" of the peaks
        s = abs(diff(ps(:,2)));
        % Use standard deviation like [Chen2011] to determine if we're in 1
        % or 2.
        x = std(s/max(s));
        
        threshold = 0.1;
        if x < threshold
            c = 1;
        else
            c = 2;
        end 
    end
end
function ps = get_turning_points(y)
    diff_s = sign(diff(y));
    % diffs(1) = sign(x(2) - x(1))
    last_diff_sign = 0;
    ps = [];
    for i = 1:length(diff_s)
        if last_diff_sign == 0
            % Set last_diff_sign = diff_s(i) until it's nonzero
            last_diff_sign = diff_s(i);
        elseif diff_s(i) == 0
            % Ignore 0-sign diffs
        elseif diff_s(i) ~= last_diff_sign
            % last_diff_sign and diff_s are nonzero, and differ
            % We've hit a turning point!
            ps(end+1,:) = [i, y(i)];
            last_diff_sign = diff_s(i);
        end
    end
end
function ls=lengths_of_one_runs(dat)
    current_run_length = 0;
    ls = [];
    for i = 1:length(dat)
        if dat(i) == 0 && current_run_length > 0
            ls(end+1) = current_run_length;
            current_run_length = 0;
        else
            current_run_length = current_run_length + 1;
        end
    end
    if current_run_length > 0
        ls(end+1) = current_run_length;
    end
end
function p=roundoff_gauss_p(expected_err)
    if abs(expected_err) > 6
        p = 0;
    else
        norm_constant = 1;
        sigma_2 = 0.8;
        p = norm_constant * exp(-expected_err*expected_err/(2*sigma_2));
    end
end
function bs=blocks_8x8(dat)
    s = size(dat);
    n_blocks = s/8;
    % Create a 1xn_blocks(1) vector of 8s
    blocks_x = repmat(8, 1, n_blocks(1));
    blocks_y = repmat(8, 1, n_blocks(2));
    % Get a 2D array of blocks
    block_2d_array = mat2cell(dat, blocks_x, blocks_y);
    % Convert to a linear array of blocks
    bs = reshape(block_2d_array, [1 numel(block_2d_array)]);
end
% Functions for importing a png with [0,1] range color values
function dat=import_png_bw(name)
    I = imread(name);
    dat = (I);
end
function [Y,Cr,Cb]=import_png_ycrcb(name)
    I = imread(name);
    % Get combined full-res YCrCb
    YCrCb = rgb2ycbcr((I));
    % Split off Y component
    Y = YCrCb(:,:,1);
    % Resample Cr, Cb at 1/2 resolution
    Cr = YCrCb(1:2:end,1:2:end,2);
    Cb = YCrCb(1:2:end,1:2:end,3);
end

% From https://github.com/rygorous/dct_blog, https://fgiesen.wordpress.com/2013/11/04/bink-2-2-integer-dct-design-part-1/
% function out=bink_idct_B2_partial(in,stages)
%   % extract rows (with input permutation)
%   c0 = in(1,:);
%   d4 = in(2,:);
%   c2 = in(3,:);
%   d6 = in(4,:);
%   c1 = in(5,:);
%   d5 = in(6,:);
%   c3 = in(7,:);
%   d7 = in(8,:);
% 
%   % odd stage 4
%   c4 = d4;
%   c5 = d5 + d6;
%   c7 = d5 - d6;
%   c6 = d7;
% 
%   if stages > 1
%     % odd stage 3
%     b4 = c4 + c5;
%     b5 = c4 - c5;
%     b6 = c6 + c7;
%     b7 = c6 - c7;
% 
%     % even stage 3
%     b0 = c0 + c1;
%     b1 = c0 - c1;
%     b2 = c2 + c2/4 + c3/2;
%     b3 = c2/2 - c3 - c3/4;
% 
%     if stages > 2
%       % odd stage 2
%       a4 = b7/4 + b4 + b4/4 - b4/16;
%       a7 = b4/4 - b7 - b7/4 + b7/16;
%       a5 = b5 - b6 + b6/4 + b6/16;
%       a6 = b6 + b5 - b5/4 - b5/16;
% 
%       % even stage 2
%       a0 = b0 + b2;
%       a1 = b1 + b3;
%       a2 = b1 - b3;
%       a3 = b0 - b2;
% 
%       if stages > 3
%         % stage 1
%         o0 = a0 + a4;
%         o1 = a1 + a5;
%         o2 = a2 + a6;
%         o3 = a3 + a7;
%         o4 = a3 - a7;
%         o5 = a2 - a6;
%         o6 = a1 - a5;
%         o7 = a0 - a4;
%       else
%         o0 = b0;
%         o1 = b1;
%         o2 = b2;
%         o3 = b3;
%         o4 = b4;
%         o5 = b5;
%         o6 = b6;
%         o7 = b7;
%       end
%     else
%       o0 = b0;
%       o1 = b1;
%       o2 = b2;
%       o3 = b3;
%       o4 = b4;
%       o5 = b5;
%       o6 = b6;
%       o7 = b7;
%     end
%   else
%     o0 = c0;
%     o1 = c1;
%     o2 = c2;
%     o3 = c3;
%     o4 = c4;
%     o5 = c5;
%     o6 = c6;
%     o7 = c7;
%   end
% 
%   % output
%   out = [o0; o1; o2; o3; o4; o5; o6; o7];
% end