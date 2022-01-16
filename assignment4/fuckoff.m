Di = image_DCT_params("test1.png", "gray");
x = abs(Di{1,2});
h = histcounts(abs(x),'BinMethod','Integer');

x_q_5 = 5 * floor(x ./ 5);
h_q_5 = histcounts(abs(x_q_5),'BinMethod','Integer');
h_q_5_prime = filter([1 1 1 1 1], 1, h);
% h_q_5_prime = h_q_5_prime(2:end);
h_q_5_prime(2:5:end) = 0;
h_q_5_prime(3:5:end) = 0;
h_q_5_prime(4:5:end) = 0;
h_q_5_prime(5:5:end) = 0;
figure;
plot(h_q_5);
hold on;
plot(h_q_5_prime);
hold off;
title("Comparing h_q_5 computed from quantization to filter+sample");

x_q_5_q_11 = 11 * floor(x_q_5 ./ 11);
h_q_5_q_11 = histcounts(abs(x_q_5_q_11),'BinMethod','Integer');

x_q_5_q_3 = 3 * floor(x_q_5 ./ 3);
h_q_5_q_3 = histcounts(abs(x_q_5_q_3),'BinMethod','Integer');


x_q_5_q_10 = 10 * floor(x_q_5 ./ 10);
h_q_5_q_10 = histcounts(abs(x_q_5_q_10),'BinMethod','Integer');
% h_q_5_q_10 = h_q_5_q_10(1:length(h));

x_q_10 = 10 * floor(x ./ 10);
h_q_10 = histcounts(abs(x_q_10),'BinMethod','Integer');

figure;
y = abs(fft(h));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x)");

figure;
y = abs(fft(filter([1 1 1 1 1], 1, h)));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(filtered h)");

figure;
y = abs(fft(h_q_5));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x quantized 5)");

figure;
y = abs(fft(h_q_5_q_3));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x quantized 5 quantized 3)");

figure;
y = abs(fft(h_q_5_q_11));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x quantized 5 quantized 11)");

figure;
y = abs(fft(h_q_5_q_10));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x quantized 5 quantized 10)");

figure;
y = abs(fft(h_q_10));
y = fftshift(y);
centre = floor(length(y)/2);
plot(centre-(0:(length(y)-1)), y);
title("fft(x quantized 10)");