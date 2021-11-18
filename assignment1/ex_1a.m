%% DSP Assignment 1 
% Samuel Stark <sws35>
%
% November 2021

%% 1.a
x = [ 0 0 0 -4 0 0 0 0 0 0 2 2 2 2 ...
2 0 -3 -3 -3 0 0 0 0 0 1 -4 0 4 ...
3 -1 2 -3 -1 0 2 -4 -2 1 0 0 0 3 ...
-3 3 -3 3 -3 3 -3 3 -3 0 0 0 0 0 0 ];
n = 0:length(x)-1;

%% 1.a) 4 Point Moving Average
a=[1];
b=[1 1 1 1]/4;
y = filter(b, a, x);
figure;
plot(n, x, 'bx-', n, y, 'ro-');
ylim([-7,7]);
title("1.a) 4-point moving average");


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
