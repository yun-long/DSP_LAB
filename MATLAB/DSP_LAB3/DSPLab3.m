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
N1 = [16 32];
% generate rectangular bartlett hamming windows
RectWin12 = rectwin(N1(1));
RectWin36 = rectwin(N1(2));
BartWin12 = bartlett(N1(1));
BartWin36 = bartlett(N1(2));
HammWin12 = hamming(N1(1));
HammWin36 = hamming(N1(2));
% calculate the DTFT of the windows
[dtftRect12, wRect12] = dtft(RectWin12);
[dtftRect36, wRect36] = dtft(RectWin36);
[dtftBart12, wBart12] = dtft(BartWin12);
[dtftBart36, wBart36] = dtft(BartWin36);
[dtftHamm12, wHamm12] = dtft(HammWin12);
[dtftHamm36, wHamm36] = dtft(HammWin36);
% plot(RectWin12);
figure
subplot(3,1,1);
plot(wRect12,abs(dtftRect12),'r',wRect36,abs(dtftRect36));
subplot(3,1,2);
plot(wBart12,abs(dtftBart12),'r',wBart36,abs(dtftBart36));
subplot(3,1,3);
plot(wHamm12,abs(dtftHamm12),'r',wHamm36,abs(dtftHamm36));
% generate a sinusoids function
N2 = N1(2);
n = 0:N2-1;
X = sin(2*pi*n/N2);
% modulate it with different windows
rectX = X.*RectWin36';
bartX = X.*BartWin36';
hammX = X.*HammWin36';
% calculate the dtft of each modulated sinusoids functions
[dtftRX, wRX] = dtft(rectX);
[dtftBX, wBX] = dtft(bartX);
[dtftHX, wHX] = dtft(hammX);
% plot the modulated function and also its dtft
figure
plot(n,rectX,'r--',n,bartX,'g-',n,hammX,'b');
figure
plot(wRX,abs(dtftRX),'r--',wBX,abs(dtftBX),'g-',wHX,abs(dtftHX),'b');
% generate two sinusoids at different frequencies
n1 = 0:N1(1)-1;
n2 = 0:N1(2)-1;
% define the parameters
a  = [0.2, 0.3];
w0 = 0.2*pi;
wa1n1 = w0 + 2*pi*a(1)/N1(1);
wa2n1 = w0 + 2*pi*a(2)/N1(1);
wa1n2 = w0 + 2*pi*a(1)/N1(2);
% 
Xa1n1 = cos(w0*n1) + cos(wa1n1*n1);
Xa2n1 = cos(w0*n1) + cos(wa2n1*n1);
Xa1n2 = cos(w0*n2) + cos(wa1n2*n2);
% 
[dtftXa1n1, wXa1n1 ]= dtft(Xa1n1);
[dtftXa2n1, wXa2n1 ]= dtft(Xa2n1);
[dtftXa1n2, wXa1n2 ]= dtft(Xa1n2);
%
figure
plot(wXa1n1,abs(dtftXa1n1),'r',wXa2n1,abs(dtftXa2n1));
figure
plot(wXa1n1,abs(dtftXa1n1),'r',wXa1n2,abs(dtftXa1n2));

%% Problem 6.3 DTMF Signal Decoding
% define parameters
fs = 8000;
N  = 256;
% generate the signal for '7'
x7 = dtmfdial('7',fs);
% calculate the dft of x7
[dftX7, wX7] = dft(x7);
yx7 = abs(dftX7);
f = wX7/(2*pi)*fs;
% plot the signal
plot(f,yx7);
hold on
fx = find(yx7 > max(yx7)*9/10);
plot(f(fx(1:2)),yx7(fx(1:2)),'ro');
%
load mynumber.mat
p = 1;
for i = 1:11
    key(i) = dtmfcoef(xx(p:p+255),fs);
    p = p + 320;
end

%% Problem 6.4 Inverse Filtering for Speech Dereverberation
% generate a room impulse response
fs  = 44100;
mic = [19, 18, 1.6];
n   = 12;
r   = 0.8;
rm  = [10, 14, 20];
src = [5, 2, 1];
h = rir(fs,mic,n,r,rm,src);
% import the speech signal 
[x, fs] = audioread('speech.wav');
%soundsc(x);
% filter the speech using generated room impulse response
y = conv(h,x);
% length(y) == length(h) + length(x) - 1;
% soundsc(y);
yx = fft(x);
nx = 0:1/length(yx):1-1/length(yx);
yy = fft(y);
ny = 0:1/length(yy):1-1/length(yy);
figure
plot(nx,abs(yx),ny,abs(yy),'r.');





