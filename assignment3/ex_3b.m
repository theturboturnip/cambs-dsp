%% 3.b.i
% Load data

f = fopen('iq-fm-96M-240k.dat', 'r', 'ieee-le');
c = fread(f, [2,inf], '*float32');
fclose(f);
z = c(1,:) + 1j*c(2,:);

% Fs = 240kHz, but signal itself was originall 96MHz +- (192/2)kHz
Fs = 240000;

%%
% TODO FM demodulate


%%
% Low-pass filter
butter()
