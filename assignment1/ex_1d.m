%% 1.d)
% Simulate the reconstruction of a sampled band-pass signal

%%
% Generate a 1 s noise sequence r, as in part (c)(i ), but this time use
% a sampling frequency of 3 kHz. Set the first and last 500 samples to
% zero.
f_s = 3000;
t = 1; % 1 second long
N = f_s * t;
ts = linspace(1/f_s, t, N);
r = randn(N,1);
r(1:500) = 0;
r(end-500:end) = 0;

%%
% Apply a band-pass filter that attenuates frequencies outside 31-44Hz
f1 = 31;
f2 = 44;
% 3 = 3rd order filter
% 30 = -30dB for frequencies outside the range
[b, a] = cheby2(3, 30, [f1 f2]/(f_s/2));
% Apply the filter
x = filtfilt(b, a, r);

%%
% Plot r, x together. r is very high frequency, so the band-pass will
% remove a lot and x will have a much lower amplitiude
figure;
hold on;
plot(r);
plot(x);
legend("r","x");
hold off;

%%
% Sample x at 30Hz, set all but every 100th value to 0
y = x;
for i = 1:99
    y(i:100:end) = 0;
end

figure;
plot(x);
hold on;
stem(y);
hold off;

%%
% Reconstruct y with sinc interpolation (see 1.c.i for working)
% Change the scaling factor to 1/100 instead of 1/3, because we sampled
% every 100th value
z = zeros(N,1);
figure;
stem(y)
hold on;
for i_y = 1:N
    data = y(i_y) * sinc((ts/t - i_y/N) * f_s/100);
    plot(data);
    z = z + data';
end
hold off;
title("1.d Per-sample sinc functions")

%%
% Generate another band-pass filter for 30-45Hz, apply to y to reconstruct
% as u. Multiply by 100 to compensate for energy lost during sampling.
f1 = 30;
f2 = 45;
% 3 = 3rd order filter
% 30 = -30dB for frequencies outside the range
[b, a] = cheby2(3, 30, [f1 f2]/(f_s/2));
u = 100 * filtfilt(b, a, y);

%% 1.d.i) Comparisons
% Plot x, y, z, and u on top of each other in one figure
figure;
plot(x);
hold on;
stem(y);
plot(z);
plot(u);
legend("x", "y", "z", "u");
title("1.d.i) x, y, z, and u")
hold off;

%% 
% x and z are not equal - sinc interpolation can only reconstruct
% frequencies up to $f_{sample}/2 = 30/2 = 15Hz$, but x contains
% higher-frequency data.
figure;
plot(x);
hold on;
stem(y);
plot(z);

legend("x", "y", "z");
xlim([1000 1500]);
title("1.c.i) x vs. z (sinc-interpolated sampled x)");
hold off;

%% 
% x and u are similar, but not equal.
% TODO - is this to do with aliasing introduced by sampling?
figure;
plot(x);
hold on;
stem(y);
plot(u);

legend("x", "y", "u");
xlim([1000 1500]);
title("1.c.i) x vs. u (low-pass filtered sampled x)");
hold off;


%% 1.d.ii)
% Q: Why does the reconstructed waveform differ much more from the original
% if you reduce the cut-off frequencies of all band-pass filters by 5 Hz?
%
% TODO