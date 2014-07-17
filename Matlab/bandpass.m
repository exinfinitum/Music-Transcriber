clear
clc
%h(t) = 2bsinc(2bt). Frequency response is rect(f/2b).
N = 10000% Number of samples in input signal.

low = -10; %Lower bound of input signal
hi = 10; %Upper bound of input signal

lowf = 10; %Limit frequency for low-pass. Freqs above this are rejected.
highf = 1; %Limit freq for high-pass. Freqs below this are rejected.

klow = -10; %Lower bound of the sinc kernel
khi = 10; %Upper bound of the sinc kernel
ksize = 10000 %Number of samples (discrete data points) in the sinc kernel

%As a rule of thumb, the kernel should be approximately equal in size
%and sample count as the piece of the input signal that we deal with (the
%buffer).

space = linspace(low, hi, N);
g = zeros(ksize, 1);
f = zeros(N, 1);
kspace = linspace(klow, khi, ksize);

for t = 1:ksize

g(t) = lowf*sinc(lowf*kspace(t)) - highf*sinc(highf*kspace(t));

end
g = g .* blackman(ksize);%Reduces the annoying errors near the tails of the result
for t = 1:N;
%g(t) = sinc(b*space(t));


%convolve f, a sin function, with g.
f(t) = sin(12.2*pi*space(t)) + sin(2*pi*space(t)) + sin(0.5*pi*space(t));% + sin(2425*pi*space(t)) ;%Fourier combination of 3 sin functions is input
%f(t) = dirac(space(t));
end

%Let h be the convolution of g and f.
h = conv(f, g, 'same');
figure
plot(h)
title('Output')
figure
plot(f)
title('Input')
