%% 3.b.i
% Load data
f = fopen('iq-fm-96M-240k.dat', 'r', 'ieee-le');
c = fread(f, [2,inf], '*float32');
fclose(f);
z = c(1,:) + 1j*c(2,:);

% Fs = 240kHz, but signal itself was originally 96MHz +- (192/2)kHz
Fs = 240000;

%%
% FM demodulate (see function)
s = fm_demodulate(Fs, z);
plot(s)

%%
% Low-pass filter
F_low = 16000;
[b_16,a_16] = butter(7, (F_low/(Fs/2)), 'low');
h = fvtool(b_16, a_16);
s_filtered = filter(b_16, a_16, s);

%%
% Reduce from 240kHz to 48kHz
s_reduced = s_filtered(1:5:end);

%%
% Normalize, write out to audio
s_normal = normalize(s_reduced, 'range', [-1 1]);
audiowrite('ex_3b_i.wav', s_normal, 48000);

%% 3.b.ii
% Load data
f = fopen('iq-fm-97M-3.6M.dat', 'r', 'ieee-le');
c = fread(f, [2,inf], '*float32');
fclose(f);
z = c(1,:) + 1j*c(2,:);

% Fs = 3.6MHz, but signal itself was originally 97MHz +- (2.88/2)MHz
Fs = 3.6E6;

Fbase = 97E6;
F1 = 96E6;
F2 = 97.2E6;
F3 = 98.5E6;

%%
% Shift the frequency up by 1MHz
z_shifted = shift_freq(z, F1 - Fbase, Fs);


%% 
% Apply a low-pass filter. 
F_low = 200E3;
[b_200,a_200] = butter(5, (F_low/(Fs/2)), 'low');
figure;
h = fvtool(b_200, a_200);
z_filtered = filter(b_200, a_200, z_shifted);

% Also, prepare a 19kHz notch filter - each signal has a 19kHz spike
[b_n, a_n] = iirnotch(19000/24000, (19/24)/35);

% Reset the 16kHz low-pass filter to be useful under 3.6MHz
[b_16,a_16] = butter(3, (16000/(Fs/2)), 'low');


%%
% Display spectrograms
%
% TODO How does the displayed frequency relate to the radio frequency
figure;
spectrogram(z,'yaxis');
title("original signal, centred on 97MHz");
% 96MHz < 97MHz, so it's in the "negative frequency" range and have to wrap
% around
yline((2+(F1 - Fbase)/(Fs/2)), '', "96MHz (\Delta = -1MHz)",'FontWeight', 'bold');
% F2 and F3 > 97 => they're in the positive frequency range
yline(([F2 F3] - Fbase)/(Fs/2), '', ["97.2MHz (\Delta = 0.2MHz)", "98.5MHz (\Delta = 1.5MHz)"],'FontWeight', 'bold');
% Show the split between +ve, -ve frequencies
yline(1, '--', 'Negative-Positive Frequency Cutoff', 'LabelHorizontalAlignment', 'left', 'LineWidth', 2);


figure;
spectrogram(z_shifted,'yaxis');
title("signal shifted by 1MHz");
figure;
spectrogram(z_filtered,'yaxis');
title("signal shifted by 1MHz, low-pass filtered");

%%
% demodulate, low-pass filter, subsample, output as wav
s = fm_demodulate(Fs, z_filtered);
% Reduce from 3.6MHz to 48kHz (every 75th)
s_filtered = filter(b_16, a_16, s);
% Remove 19kHz spike
s_reduced = filter(b_n, a_n, s_reduced);
s_reduced = s_filtered(1:75:end);
s_normal = normalize(s_reduced, 'range', [-1 1]);
% s_normal = filter(b_n, a_n, s_normal);
% fft_plot(s_normal, 48000);
audiowrite('ex_3b_ii_1.wav', s_normal, 48000);

%%
% demodulate etc. F2

% Shift base frequency to F2
z_shifted = shift_freq(z, (F2 - Fbase), Fs);
figure;
spectrogram(z_shifted,'yaxis');
title("z shifted 2");
z_filtered = filter(b_200, a_200, z_shifted);
figure;
spectrogram(z_filtered,'yaxis');
title("z filtered 2");
s = fm_demodulate(Fs, z_filtered);
figure;
spectrogram(s,'yaxis');
title("s 2");
s_filtered = filter(b_16, a_16, s);
% Reduce from 3.6MHz to 48kHz (every 75th)
s_reduced = s_filtered(1:75:end);
% Remove 19kHz spike
s_reduced = filter(b_n, a_n, s_reduced);
s_normal = normalize(s_reduced, 'range', [-1 1]);
% s_normal = filter(b_n, a_n, s_normal);
% fft_plot(s_normal, 48000);
audiowrite('ex_3b_ii_2.wav', s_normal, 48000);

%%
% demodulate etc. F3

% Shift base frequency to F3
z_shifted = shift_freq(z, (F3 - Fbase), Fs);
z_filtered = filter(b_200, a_200, z_shifted);
s = fm_demodulate(Fs, z_filtered);
s_filtered = filter(b_16, a_16, s);
% Reduce from 3.6MHz to 48kHz (every 75th)
s_reduced = s_filtered(1:75:end);
% Remove 19kHz spike
s_reduced = filter(b_n, a_n, s_reduced);
s_normal = normalize(s_reduced, 'range', [-1 1]);
fft_plot(s_normal, 48000);
% Remove 19kHz signal
% s_normal = filter(b_n, a_n, s_normal);
% fft_plot(s_normal, 48000);
audiowrite('ex_3b_ii_3.wav', s_normal, 48000);

%% Functions

function s = fm_demodulate(Fs, z)
    dt = 1/Fs;
    s = angle(z(2:end)./z(1:(end-1)))/dt;

%     dt = 1/Fs;
%     dz = z(2:end) - z(1:(end-1));
%     z_short = z(1:(end-1));
%     s = imag((dz/dt).*conj(z_short))./(abs(z_short).^2);
%     s = ((dz/dt)./z(2:end))/(2i*pi*1);
end

function z_prime=shift_freq(z, f, Fs)
    phasor = exp(2*pi*-1i*f*(1:length(z))/Fs);
    z_prime = z .* phasor;
end

function fft_plot(x, Fs)
    y = fft(x);
    % Assume symmetry
    y_single = abs(y(1:(size(y,2)/2)))*2;
    figure;
    plot(linspace(0,Fs/2,size(y_single,2)), y_single);
end