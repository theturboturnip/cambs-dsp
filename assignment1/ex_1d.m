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

%% 1.d.i) x vs z
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
title("1.d.i) x vs. z (sinc-interpolated sampled x)");
hold off;

%% 1.d.i) x vs u
% x and u are similar, but not equal.
% I think this is because the band-pass filter is not perfect.
figure;
plot(x);
hold on;
stem(y);
plot(u);

legend("x", "y", "u");
xlim([1000 1500]);
title("1.d.i) x vs. u (30-45Hz band-pass filtered sampled x)");
hold off;

%%
% Frequency response of the second band-pass filter.
% Note that the values for the 31-44Hz range (used to generate x) drop to
% around -20dB. These frequencies are unfairly reduced when generating u,
% so it doesn't match x exactly.
f1 = 30;
f2 = 45;
[b, a] = cheby2(3, 30, [f1 f2]/(f_s/2));
figure;
freqz(b,a,2048,f_s);
xlim([20, 50]);
xline(31);
xline(44);


%% 1.d.ii)
% Q: Why does the reconstructed waveform differ much more from the original
% if you reduce the cut-off frequencies of all band-pass filters by 5 Hz?
%
% The question specifically states that the band-pass filters change, but
% does not state that the sampling frequency for y changes from 30Hz.
%
% When the band-pass frequencies change, x's frequency range is 26-39Hz.
% The sampling frequency for y, 30Hz, is within this range.
% This causes aliasing, which distorts the frequency domain representation
% (like in slide 58) and thus distorts the output of the low-pass filter.
% 
% Before, x's range was 31-44Hz, and f_y = 30Hz was outside this range, so
% there wasn't a problem.
% With the decreased frequencies we can also resolve the problem by
% decreasing f_y to 25Hz, and the final figure shows this fixes the issue.

% Find x' = x but frequencies reduced by 5
f1 = 31-5;
f2 = 44-5;
[b, a] = cheby2(3, 30, [f1 f2]/(f_s/2));
x_prime = filtfilt(b, a, r);

% Find y'_30 = sampling of x' @ 30Hz
y_prime_30 = x_prime;
for i = 1:99
    y_prime_30(i:100:end) = 0;
end
% Find y'_25 = sampling of x' @ 25Hz
y_prime_25 = x_prime;
for i = 1:119
    y_prime_25(i:120:end) = 0;
end

% Find u' = reconstruction of x' from filter of y' with freqs reduced by 5
f1 = 30-5;
f2 = 45-5;
[b, a] = cheby2(3, 30, [f1 f2]/(f_s/2));
u_prime_30 = 100 * filtfilt(b, a, y_prime_30);
u_prime_25 = 100 * filtfilt(b, a, y_prime_25);

% Plot all
figure;
plot(x_prime);
hold on;
stem(y_prime_30);
plot(u_prime_30);
legend("x'", "y' @ 30Hz", "u' from y' @ 30Hz");
title("1.d.ii) x' vs. u'\\(band-pass frequencies reduced by 5Hz, sampled @ 30Hz)");
hold off;

figure;
plot(x_prime);
hold on;
stem(y_prime_25);
plot(u_prime_25);
legend("x'", "y' @ 25Hz", "u' from y' @ 25Hz");
title("1.d.ii) x' vs. u'\\(band-pass frequencies reduced by 5Hz, sampled @ 25Hz)");
hold off;
