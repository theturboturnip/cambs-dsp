%% 3.b.i
% Load data

f = fopen('iq-fm-96M-240k.dat', 'r', 'ieee-le');
c = fread(f, [2,inf], '*float32');
fclose(f);
z = c(1,:) + 1j*c(2,:);

% Fs = 240kHz, but signal itself was originall 96MHz +- (192/2)kHz
Fs = 240000;

%%
% FM demodulate (see function)
s = fm_demodulate(Fs, z);
plot(s)

%%
% Low-pass filter
F_low = 16000;
[b,a] = butter(10, (F_low/(Fs/2)), 'low');
h = fvtool(b, a);
s_filtered = filter(b, a, s);

%%
% Reduce from 240kHz to 48kHz
s_reduced = s_filtered(1:5:end);

%%
% Normalize, write out to audio
s_normal = normalize(s_reduced, 'range', [-1 1]);
audiowrite('ex_3b_i.wav', s_normal, 48000);

%% Functions

function s = fm_demodulate(Fs, z)
    dt = 1/Fs;
    s = angle(z(2:end)./z(1:(end-1)))/dt;
end

function z_prime=shift_freq(z, f, fs)

end