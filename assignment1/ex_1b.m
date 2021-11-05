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
