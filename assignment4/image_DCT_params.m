function Di=image_DCT_params(name,color_type)
    if color_type == "color"
        [Y,Cr,Cb] = import_png_ycrcb(name);
        dat = Y;
    else
        dat = import_png_bw(name);
    end
    
    bs = blocks_8x8(dat);
    Di = DCT_params(bs);
end

function Di=DCT_params(bs)
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
%         if ismember(255, b) || ismember(0, b) || length(unique(b)) == 1
%             % This block is potentially-truncated, or uniform
%             % => DCT could throw off results
%         else
            for j = 1:8
                for k = 1:8
                    % TODO
                    Di{j,k}(i) = d(j,k);
                end
            end
%         end
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