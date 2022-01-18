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
    % Di = 8x8 cell array, each cell = 1xlength(bs) vector
    Di = cell(8, 8);
    % Preallocate enough memory for all blocks on each coefficient
    for j = 1:8
        for k = 1:8
            Di{j,k} = zeros(1, length(bs));
        end
    end
    % Count the number of DCTs we actually compute
    n_dcts = 1;
    for i = 1:length(bs)
        b = bs{i};
        
        % Do Fan2003 filtering here
        if ismember(255, b) || ismember(0, b) || length(unique(b)) == 1
            % This block is potentially-truncated, or uniform
            % => DCT could throw off results, don't do it
        else
            % Perform a DCT
            % For unsigned 8-bit input, DCT coefficients are signed 11-bit ints
            % => it's OK to round them
            d = round(dct2(b));
            for j = 1:8
                for k = 1:8
                    Di{j,k}(n_dcts) = d(j,k);
                end
            end
            n_dcts = n_dcts + 1;
        end
    end
    % Shrink the individual cell vectors to only the DCTs we calculated
    for j = 1:8
        for k = 1:8
            Di{j,k} = Di{j,k}(1:(n_dcts-1));
        end
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