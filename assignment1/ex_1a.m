%% 1.a
x = [ 0 0 0 -4 0 0 0 0 0 0 2 2 2 2 ...
2 0 -3 -3 -3 0 0 0 0 0 1 -4 0 4 ...
3 -1 2 -3 -1 0 2 -4 -2 1 0 0 0 3 ...
-3 3 -3 3 -3 3 -3 3 -3 0 0 0 0 0 0 ];
n = 0:length(x)-1;

a=[1];
b=[1 1 1 1]/4;
y = filter(b, a, x);
figure;
plot(n, x, 'bx-', n, y, 'ro-');
ylim([-7,7]);
title("Example - 4-point moving average");


%% 1.a.i) Exponential Averaging System
% $y_n = alpha * x_n + (1 - alpha)*y_{n-1}$
%
% $=> y_n - (1 - alpha) * y_{n-1} = alpha * x_n$
%
% $=> [1, -(1 - alpha)] . [y_n, y_{n-1}] = [alpha] . [x_n]$
alpha = 0.5;
a=[1, -(1 - alpha)];
b=[alpha];
y = filter(b, a, x);
figure;
plot(n, x, 'bx-', n, y, 'ro-');
ylim([-7,7]);
title("1.a.i) Exponential Averaging System ({\alpha = 0.5})");

%% 1.a.ii) Accumulator
% $$y_n = x_n + y_{n-1}$$
%
% $$=> y_n - y_{n-1} = x_n$$
%
% $$=> [1, -1] . [y_n, y_{n-1}] = [1] . [x_n]$$
a=[1, -1];
b=[1];
y = filter(b, a, x);
figure;
plot(n, x, 'bx-', n, y, 'ro-');
ylim([-7,7]);
title("1.a.ii) Accumulator");

%% 1.a.iii) Backwards Difference System
% $$y_n = x_n - x_{n-1}$$
%
% $$[1] . [y_n] = [1, -1] . [x_n, x_{n-1}]$$
a=[1];
b=[1, -1];
y = filter(b, a, x);
figure;
plot(n, x, 'bx-', n, y, 'ro-');
ylim([-7,7]);
title("1.a.iii) Backwards Difference");

%% 1.b
% Get table date, taking only the 'date' and 'newCasesByPublishDate'
% columns
newcases = readtable("covid_new_cases_2021_11_04.csv", 'Range', 'D:E');
% Sort table by 'date'
newcases = sortrows(newcases, 1);
% Take x = newcases.newCasesByPublishDate
x = newcases.newCasesByPublishDate;

%% 1.b - 7 Point Moving Average
a=[1];
b=[1 1 1 1 1 1 1]/7;
y = filter(b, a, x);
figure;
plot(newcases.date, y)
title("1.b Covid Data - 7-point moving average");

%% 1.b.i) Exponential Averaging System
% $y_n = alpha * x_n + (1 - alpha)*y_{n-1}$
%
% $=> y_n - (1 - alpha) * y_{n-1} = alpha * x_n$
%
% $=> [1, -(1 - alpha)] . [y_n, y_{n-1}] = [alpha] . [x_n]$
alpha = 0.5;
a=[1, -(1 - alpha)];
b=[alpha];
y = filter(b, a, x);
figure;
plot(newcases.date, y)
title("1.b.i) Exponential Averaging System ({\alpha = 0.5})");

%% 1.b.ii) Accumulator
% $$y_n = x_n + y_{n-1}$$
%
% $$=> y_n - y_{n-1} = x_n$$
%
% $$=> [1, -1] . [y_n, y_{n-1}] = [1] . [x_n]$$
a=[1, -1];
b=[1];
y = filter(b, a, x);
figure;
plot(newcases.date, y)
title("1.b.ii) Accumulator");

%% 1.b.iii) Backwards Difference System
% $$y_n = x_n - x_{n-1}$$
%
% $$[1] . [y_n] = [1, -1] . [x_n, x_{n-1}]$$
a=[1];
b=[1, -1];
y = filter(b, a, x);
figure;
plot(newcases.date, y)
title("1.b.iii) Backwards Difference");


%% 1.c
% Generate a one second long Gaussian noise sequence r with a sampling rate
% of 300Hz
f_s = 300;
t = 1; % 1 second long
N = f_s * t;
ts = linspace(0, t, N);
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
% Sample x at 100Hz by setting all but every third sample value to zero
y = x;
% element 1, 4, 7... = 0
y(1:3:end) = 0;
% element 2, 5, 8... = 0
y(2:3:end) = 0;
% element 3, 6, 9... = unchanged

%%
% Implement sinc interpolation to reconstruct the zeroed samples of y
%
% $$x(t) = \sum_{n=-\infty}^{\infty} x_n * sinc(t/t_s - n)$$
%
% MATLAB version:
% $$z = \sum_{i_y=1}^{N} y(i_y) * sinc((ts/t - i_y/N) * f_s/3)$$
%
% Translate the sinc by i_y/N to align it with the sample,
% and scale it by f_s over 3 to align the zero crossings with the other
% sample points
z = zeros(N,1);
figure;
stem(y)
hold on;
for i_y = 1:N
    data = y(i_y) * sinc((ts/t - i_y/N) * f_s/3);
    plot(data);
    z = z + data';
end
hold off;

%%
% Generate another low-pass filter with fir1, cut-off frequency 50Hz
f_c = 50; % cutoff frequency
a = [1]; % FIR has no y-terms
b = fir1(50, f_c/(f_s/2));
% Apply it to y, resulting in interpolated sequence u.
% multiply by 3 to compensate for energy lost during sampling.
u = 3 * filtfilt(b, a, y);


%%
% Plot x, y, z, and u on top of each other in one figure 
plot(x);
hold on;
stem(y);
plot(z);
plot(u);
hold off;

%% 
% TODO Compare x, y, z, u