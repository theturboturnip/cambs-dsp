%% 1.c
% Generate a one second long Gaussian noise sequence r with a sampling rate
% of 300Hz
f_s = 300;
t = 1; % 1 second long
N = f_s * t;
ts = linspace(1/f_s, t, N);
r = randn(N,1);
% Taper $r$ by setting its first and last 15 samples to zero.
r(1:15) = 0;
r(end-15:end) = 0;

%%
% Make a Finite-Impulse Response low-pass filter with cut-off frequency
f_c = 45; % cutoff frequency
a = [1]; % FIR has no y-terms
b = fir1(50, f_c/(f_s/2));
% Use the filtfilt function in order to apply that filter to the
% generated noise signal, resulting in the filtered noise signal x
x = filtfilt(b, a, r);

%%
% Plot r, x together.
figure;
hold on;
plot(r);
plot(x);
legend("r","x");
title("1.c) Raw noise $r$ vs band-pass filtered $x$");
hold off;

%%
% Sample x at 100Hz by setting all but every third sample value to zero
y = x;
% element 1, 4, 7... = 0
y(1:3:end) = 0;
% element 2, 5, 8... = 0
y(3:3:end) = 0;
% element 3, 6, 9... = unchanged

figure;
plot(x);
hold on;
stem(y);
title("1.c) x sampled at 30Hz");
hold off;

%%
% Implement sinc interpolation to reconstruct the zeroed samples of y
%
% $$x(t) = \sum_{n=-\infty}^{\infty} x_n * sinc(t/t_s - n)$$
%
% MATLAB version:
%
% $$z = \sum_{i_y=1}^{N} y(i_y) * sinc((ts/t - i_y/N) * f_s/3)$$
%
% Translate the sinc by $i_y/N$ to align it with the sample,
% and scale it by $f_s$ over 3 to align the zero crossings with the other
% sample points.
z = zeros(N,1);
figure;
hold on;
for i_y = 1:N
    data = y(i_y) * sinc((ts/t - (i_y)/N)* f_s/3);
    plot(data);
    z = z + data';
end
stem(y)
hold off;
title("1.c) Per-sample sinc functions")

%%
% Generate another low-pass filter with fir1, cut-off frequency 50Hz
f_c = 50; % cutoff frequency
a = [1]; % FIR has no y-terms
b = fir1(50, f_c/(f_s/2));
% Apply it to y, resulting in interpolated sequence u.
% multiply by 3 to compensate for energy lost during sampling.
u = 3 * filtfilt(b, a, y);


%% 1.c.i) Comparisons
% Plot x, y, z, and u on top of each other in one figure 
figure;
plot(x);
hold on;
stem(y);
plot(z);
plot(u);
legend("x", "y", "z", "u");
title("1.c.i) x, y, z, and u")
hold off;

%% 1.c.i) x vs z
% x and z are equal
figure;
plot(x);
hold on;
stem(y);
plot(z);

legend("x", "y", "z");
xlim([100 150]);
title("1.c.i) x vs. z", "(sinc-interpolated sampled x)");
hold off;

%% 1.c.i) x vs u
% x and u are similar, but not equal.
% I think this is because the low-pass filter isn't perfect.
figure;
plot(x);
hold on;
stem(y);
plot(u);

legend("x", "y", "u");
xlim([100 150]);
title("1.c.i) x vs. u", "(low-pass filtered sampled x)");
hold off;

%%
% Frequency response of the second low-pass filter.
% The magnitude at the original cutoff (45Hz) is -3dB, so higher
% frequencies in x were be unfairly reduced when generating u.
f_c = 50;
a = [1];
b = fir1(50, f_c/(f_s/2));
figure;
freqz(b,a,2048,f_s);
xlim([0, 50]);
xline(45);
title("1.c.i) u low-pass filter frequency response");


%% 1.c.ii)
% Q: Why should the first filter have a lower cut-off frequency than the second?
%
% If the second low-pass filter had a lower cut-off frequency than the
% first filter, it would discard higher-frequency information and
% not reconstruct it.

