%% DSP_LAB_3 Frequency-Domain Analysis Using the DFT

%% Introdution
% define a sinusoidal signal
N = 32;
n = 0:N-1;
x = 0.8*sin((2*pi/N)*n);
% compute the approximated DTFT using an FFT with M = 16N frequency points
M = 16*N;
w = 2*pi*(0:1/M:1-1/M);
X = fft(x,M);
% plot the signal in time domain and its approximated DTFT
figure
subplot(2,1,1)
stem(n,x);
subplot(2,1,2);
plot(w,abs(X));
hold on
% compute the DFT of the same sinusoidal signal, the DFT is computed using
% an FFT with M = N frequency points, where N is the signal length.
Xk = fft(x);
N  = length(Xk)
wk = 2*pi*(0:1/N:1-1/N);
% plot the DFT coefficients
plot(wk,abs(Xk),'ro');

%% Problem 6.1 Computational Complexity of the DFT/FFT
clear
N = 32;
n = 0:N-1;
x = 0.8*sin(0.2*pi*n);
% compute the dft using direct computation of the dft
Y1 = dft(x);
% compute the dft via fft function
Y2 = fft(x);
% plot 
figure
subplot(2,1,1);
stem(abs(Y1));
subplot(2,1,2);
stem(abs(Y2));
% compare the different results by using different N
N = [2^4 2^5 2^6 2^7 2^8];
for k = 1:5
    n = 1:N(k)-1;
    x = 0.8*sin(0.2*pi*n);
    tic;
    y1 = dft(x);
    rd(k) = toc;
    tic
    y2 = fft(x);
    tf(k) = toc;
    subplot(5,1,k);
    stem(abs(y1));
end
%% Problem 6.2 Discrete-Time Fourier Transform

