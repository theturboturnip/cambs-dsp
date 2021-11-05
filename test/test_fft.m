pkg load signal

f_s = 160E3;
f1 = 23E3;
a1 = 1;
f2 = 25E3;
a2 = 1;
block_size = 512;
blocks_per_s = f_s / block_size;
t_per_block = block_size / f_s;
n_blocks = 5;

t = 0:(1/f_s):(t_per_block * n_blocks);
sampled_signal = a1*sin(2*pi*f1*t) + a2*sin(2*pi*f2*t);

block_0 = sampled_signal(1:block_size);
t_block_0 = t(1:block_size);
plot(t_block_0, block_0);

fft_block_0 = fft(block_0 .* hamming(block_size)');
plot(fft_block_0);

# Compute the two-sided spectrum P2.
# Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
Y=fft_block_0;
L=block_size;
P2 = abs(Y/L); # Normalize impact of FFT length
# P2 is two-sided eg. includes negative frequencies
# => one-side = (first half of P2?) * 2
P1 = P2(1:L/2+1); 
P1(2:end-1) = 2*P1(2:end-1);

# Correction step - 4d
# 25kHz will have a higher peak than 23kHz, because ??? "scalloping"
# NOT the integral of sin(25kHz)^2 will cover more area.

#w1 = 2*pi*f1
## integral of sin^2(wx) = x/2-(sin(2xw)/4w)+C
## x = 0 to t_per_block
#area1 = (t_per_block/2 - sin(2*t_per_block*w1)/(4*w1)) - 0
#w2 = 2*pi*f2
#area2 = (t_per_block/2 - sin(2*t_per_block*w2)/(4*w2)) - 0

b=[1 1 1 1]/4;
a=[1];
smoothed = filter(b, a, P1);

plot(P1)
hold on;
plot(smoothed)
legend("fft", "smoothed")
hold off;
# peak 23kHz = 0.75 peak 25kHz = 1

