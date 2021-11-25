%% DSP Assignment 3
% Samuel Stark (sws35)
%
% November 2021

%% 3.a
% Deconvolve stars-blurred
%
% Configuration elements follow:
edge_taper=32; % Pixel range of tapering
zero_thresh = 1E1; % Values with a magnitude below this are clamped to that

%% 
% Load the images in floating-point
% (the original samples are uint16 = [0,65535)
convolved = single(imread("stars-blurred.png"))/65534;
filter = single(imread("stars-psf.png"))/65534;

figure;
imagesc(convolved)
figure;
imagesc(filter)

%%
% Reduce the pixel values near the edges
if edge_taper > 0
    % Lerp between the image and zero near the edges
    [h,w] = size(convolved);
    factor = zeros(h,w);
    for i = 1:w
        for j = 1:h
            % Compute minimum distance from edge
            % If at edge, e.g. i=1,j=1, distance = 0
            min_dist_from_edge = min([
%                 (i-1)^2+(j-1)^2
%                 (w-i)^2+(j-1)^2
%                 (w-i)^2+(h-j)^2
%                 (i-1)^2+(h-j)^2
                i-1
                w-i
                j-1
                h-j
            ]);

            % If dist from edge >= edge_taper, 
            % no taper => factor = 1
            % If dist from edge = 0
            % full taper => factor = 0
            factor(j,i) = min(1, min_dist_from_edge/edge_taper);
        end
    end
    convolved_tapered = (factor .* convolved);
else
    convolved_tapered = convolved;
end
figure;
imagesc(convolved_tapered);
title("Convlved image tapered at edge");

%%
% Zero pad both images to the same size
% Should be sum(convolved image, filter) sizes padded nearest pow2
sum_size = size(convolved_tapered) + size(filter);
sum_size_padded = 2.^nextpow2(sum_size);
% Evenly zero pad on both sides
convolved_padded = zero_pad(sum_size_padded, convolved_tapered);
% Pad on one side
filter_padded = zero_pad(sum_size_padded, filter);

figure;
imagesc(convolved_padded);
title("Convolved - Padded");
figure;
imagesc(filter_padded);
title("Filter - Padded");

%%
% FFT
convolved_fft = fft2(convolved_padded);
filter_fft = fft2(filter_padded);
%convolved_fft(abs(convolved_fft) < zero_thresh) = zero_thresh;
% TODO ringing
filter_fft(abs(filter_fft) < zero_thresh) = zero_thresh;

figure;
imagesc(abs(convolved_fft));
title("Convolved - FFT");
figure;
imagesc(abs(filter_fft));
title("Filter - FFT");

%% Reconstruct
% TODO This sucks lol
deconvolved_fft = convolved_fft ./ filter_fft;
figure;
imagesc(abs(deconvolved_fft));
title("Deconvolved FFT");

deconvolved_padded = ifft2(deconvolved_fft);
% TODO round
corner = (sum_size_padded-size(convolved)-size(filter))/2;
deconvolved = deconvolved_padded(corner(1):corner(1)+size(convolved,1)-1, ...
    corner(2):corner(2)+size(convolved,2)-1);

figure;
imagesc(deconvolved);
title("Deconvolved");



%% Functions
function data = zero_pad(expected_size, i)
    % four sets of padding: top, side, side, bottom
    if any(rem(size(i), 2)) % If any elements if size(i) are odd
        data = [i zeros(size(i,1), expected_size(2)-size(i,2)); zeros(expected_size(1)-size(i,1), expected_size(2))];
    else
        TOP = zeros((expected_size(1)-size(i,1))/2, expected_size(2));
        SIDE = zeros(size(i, 1), (expected_size(2)-size(i,2))/2);
        data = [ TOP ; SIDE i SIDE ; TOP];
    end
end


