Di = image_DCT_params("test1.png", "gray");

%%
x = abs(Di{1,1});
h = histcounts(abs(x),'BinMethod','Integer');

x_q_5 = 5 * floor(x ./ 5);
h_q_5 = histcounts(abs(x_q_5),'BinMethod','Integer');
h_q_5_prime = filter([1 1 1 1 1], 1, h);
% h_q_5_prime = h_q_5_prime(2:end);
h_q_5_prime(2:5:end) = 0;
h_q_5_prime(3:5:end) = 0;
h_q_5_prime(4:5:end) = 0;
h_q_5_prime(5:5:end) = 0;
% figure;
% plot(h_q_5);
% hold on;
% plot(h_q_5_prime);
% hold off;
% title("Comparing h_q_5 computed from quantization to filter+sample");

x_q_5_q_11 = 11 * floor(x_q_5 ./ 11);
h_q_5_q_11 = histcounts(abs(x_q_5_q_11),'BinMethod','Integer');

x_q_5_q_3 = 3 * floor(x_q_5 ./ 3);
h_q_5_q_3 = histcounts(abs(x_q_5_q_3),'BinMethod','Integer');

x_q_13_q_5 = 5 * floor((13 * floor(x ./ 13)) ./ 5);
h_q_13_q_5 = histcounts(abs(x_q_13_q_5),'BinMethod','Integer');

x_q_8_q_7 = 7 * floor((8 * floor(x ./ 8)) ./ 7);
h_q_8_q_7 = histcounts(abs(x_q_8_q_7),'BinMethod','Integer');

x_q_5_q_10 = 10 * floor(x_q_5 ./ 10);
h_q_5_q_10 = histcounts(abs(x_q_5_q_10),'BinMethod','Integer');
% h_q_5_q_10 = h_q_5_q_10(1:length(h));

x_q_10 = 10 * floor(x ./ 10);
h_q_10 = histcounts(abs(x_q_10),'BinMethod','Integer');

figure;
histogram(abs(x),'BinMethod','Integer');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
xlim([300 600]);

figure;
histogram(abs(x_q_5),'BinMethod','Integer');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
xlim([300 600]);

figure('Name', 'B(f)');
y = abs(fft(h));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 2);
ylim([-2*(max(y)/2)/6, max(y)/2]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
% set(gca,'Visible','off')
% set(gca,'XColor','none')
% xline(0);
% yline(0);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;

% figure('Name', 'q*sinc(fq)');
hold on;
fs = linspace(-0.5, 0.5, length(y));
yyaxis right;
plot(fs*length(y), 5 * sinc(fs * 5), 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca, 'YColor',[0 0 0])
% xlim([-0.5 0.5]);
ylim([-2 6]);
% xline(0);
% yline(0);
hold off;
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% grid on;
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'Visible','off')
% set(gca,'XColor','none')

figure('Name', 'Dirac Comb(1/q)');
% diraccomb = zeros(size(y));
period = round(length(y)/5);
% diraccomb(centre:period:end) = 1;
% diraccomb(centre:-period:0) = 1;
diraccomb = [centre:-period:0  centre:period:length(y)];
h=stem(diraccomb - centre, 5*ones(size(diraccomb)), 'LineWidth', 2);
h.MarkerFaceColor = 'w';
% plot(centre-(0:(length(y)-1)), diraccomb);
ylim([-2 6])
xlim([0-centre, length(y)-centre])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
% set(gca,'Visible','off')
% set(gca,'XColor','none')
% xline(0);
% yline(0);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;

%% 
figure('Name', "B(f) 1x quantized (q=5)");
y = abs(fft(h_q_5));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
ylim([-2*(max(y)*1.1)/6, max(y)*1.1]);
grid on;

figure('Name', "B(f) 2x quantized (q1=5, q2=3)");
y = abs(fft(h_q_5_q_3));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;

figure('Name', "B(f) 2x quantized (q1=5, q2=11)");
y = abs(fft(h_q_5_q_11));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 2);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;

%%
x2 = abs(Di{1,2});
h2 = histcounts(abs(x2),'BinMethod','Integer');

x2_q_5 = 5 * floor(x2 ./ 5);
h2_q_5 = histcounts(abs(x2_q_5),'BinMethod','Integer');

x2_q_5_q_11 = 11 * floor(x2_q_5 ./ 11);
h2_q_5_q_11 = histcounts(abs(x2_q_5_q_11),'BinMethod','Integer');

x2_q_5_q_3 = 3 * floor(x2_q_5 ./ 3);
h2_q_5_q_3 = histcounts(abs(x2_q_5_q_3),'BinMethod','Integer');

figure('Name', 'AC B(f)');
y = abs(fft(h2));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 3);
ylim([0, max(y)*1.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;

figure('Name', "AC B(f) 1x quantized (q=5)");
y = abs(fft(h2_q_5));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 3);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
ylim([0, max(y)*1.1]);
grid on;

figure('Name', "AC B(f) 2x quantized (q1=5, q2=3)");
y = abs(fft(h2_q_5_q_3));
y = fftshift(y)./length(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y, 'LineWidth', 3);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
ax = gca;
ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
grid on;
ylim([0, max(y)*1.1]);

%%
% figure;
% y = abs(fft(h));
% y = fftshift(y)./length(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x)");
% 
% figure;
% y = abs(fft(filter([1 1 1 1 1], 1, h)));
% y = fftshift(y)./length(y);
% centre = floor(length(y)/2);
% fs = (centre-(0:(length(y)-1)) )./ length(y);
% plot(fs, y);
% hold on;
% plot(fs, 5*sinc(fs .* 5));
% hold off;
% title("fft(filtered h)");
% 
% figure;
% y = abs(fft(h_q_5));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 5)");
% 
% figure;
% y = abs(fft(h_q_5_q_13));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 5 quantized 13)");
% 
% figure;
% y = abs(fft(h_q_13_q_5));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 13 quantized 5)");
% 
% figure;
% y = abs(fft(h_q_5_q_11));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 5 quantized 11)");
% 
% figure;
% y = abs(fft(h_q_8_q_7));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 8 quantized 7)");

% figure;
% y = abs(fft(h_q_5_q_10));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 5 quantized 10)");
% 
% figure;
% y = abs(fft(h_q_10));
% y = fftshift(y);
% centre = floor(length(y)/2);
% plot(centre-(0:(length(y)-1)), y);
% title("fft(x quantized 10)");