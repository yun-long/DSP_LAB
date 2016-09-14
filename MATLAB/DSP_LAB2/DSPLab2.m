%% Diginal Singnal Processing Lab 2

%% 3.4 Filtering with Matlab
clear
n  = -10:20;
r  = 0.8;
w0 = 0.5*pi;
b  = [1  0            -1];
a  = [1  -2*r*cos(w0) r^2];
% calculate the impulse response
delta =+ (n==0);
h  = filter(b,a,delta);
stem(n,h);
% calculate the frequency response
[H,W] = freqz(b,a);
plot(W/(pi),20*log10(abs(H)));
% delta omega
deltaW = 2*(1-r)/sqrt(r);
% ? How to find the 3db, and mark it
%% Preparation
clear
% 1)
n = -10:20;
% unit step function
u0  =+ (n >= 0);
u6  =+ (n >= 6); 
% Kronecker delta function
delta1 =+ (n==1);
delta2 =+ (n==2);
delta3 =+ (n==3);
% different kind of signals
x1 = sin(0.12*pi*n);
x2 = u0 - u6;
x3 = 0.9.^n.*u0;
x4 = 0.5*delta1 + delta2+ + 0.5*delta3;
x5 = 0.9.^n.*cos(0.2*pi*n);
x6 = sinc(0.2*pi*n);
% 2) N.A.
% 3) N.A.
% 4) N.A.

% Problem 6.1 Discrete-Time Signals
figure
subplot(3,2,1);
plot(n,x1);
title('sine function');
subplot(3,2,2);
plot(n,x2);
title('unit step function');
subplot(3,2,3);
plot(n,x3);
title('exponential function');
subplot(3,2,4);
plot(n,x4);
title('Kronecker delta function');
subplot(3,2,5);
plot(n,x5);
title('cosine function');
subplot(3,2,6);
plot(n,x6);
title('sinc function');
% implement convolution function in convmat.m
close all;
n = 1:36;
w1 = 0.2*pi;
w2 = 0.4*pi;
h([1:5])   = 0.2;
x1n([1:10]) = 1;
x2n = sin(w1*n);
x3n = sin(w2*n);
y11 = convmat(h,x1n);
y12 = convmat(h,x2n);
y13 = convmat(h,x3n);
y21 = conv(h,x1n);
y22 = conv(h,x2n);
y23 = conv(h,x3n);
figure
subplot(2,1,1);
stem(y13);
subplot(2,1,2);
stem(y23);

%% Problem 6.2 Musical Tone Synthesis
n = 0:8191;
f1 = 392;
fs = 8192;
% generate a bisc signal
s1  = sin(2*pi*(f1/fs)*n);
% add 7 harmonics with differenct frequeucies to the fundamental
% oscilattion
s2 = s1;
for k = 2:8
    s2 = s2 + 0.25^(k-1)*sin(2*pi*k*(f1/fs)*n); 
end
% generate an envelope A(n)
N = length(s2);
n = [240,7200];
An = envelope(N,0.25,n);
% modulate signal s2 with the new geneated envelope A(n)
s3 = An.*s2;
% plot the signal s2 and modulated signal s3
figure
subplot(2,1,1);
plot(s2);
subplot(2,1,2);
plot(s3);
% plot it in frequency domain by plotting the FFT.
Y1 = fft(s1);
Y2 = fft(s2);
Y3 = fft(s3);
figure 
subplot(3,1,1);
plot(1:length(Y1),abs(Y1)/max(abs(Y1)));
subplot(3,1,2);
plot(1:length(Y2),abs(Y2)/max(abs(Y2)));
subplot(3,1,3);
plot(1:length(Y3),abs(Y3)/max(abs(Y3)));
% listen to the musical tone
% soundsc(s1);
% soundsc(s2);
% soundsc(s3);
% load the pianoG3.mat
load pianoG3.mat
soundsc(g);
Yg = fft(g);
figure 
% compare real recorded piano sound with synthesized tone from s3 in
% frequency domain
subplot(2,1,1);
plot(abs(Yg)/max(abs(Yg)));
subplot(2,1,2);
plot(abs(Y3)/max(abs(Y3)));

%% Problem 6.3 Discrete-Time Systems
% define the parameters for the transfer function
clear
n = -10:20;
b = [0.16  0.48  0.48  0.16];
a = [1.00  0.13  0.52  0.30];
% calculate the impulse response
delta =+ (n==0);
h = filter(b,a,delta);
% calculate the frequency response
[H,W] = freqz(b,a);
% plot the impulse and frequency respones
figure
subplot(2,1,1);
stem(h);
subplot(2,1,2);
plot(W/pi,20*log10(abs(H)));
% measure the frequency response magnitude
K = 100;
n = 0:255;
wk = pi*(0:K)./100;
% discompose the continuous frequency into discrete frequencies, in total,
% we have 100 samples, each frequency sample contains 256 time units.
for k = 1:K+1
    xk = cos((pi*k/100)*n);
    yk(k,:) = filter(b,a,xk);
    yk(k,1:30) = 0;
end
Hy = max(yk,[],2);
% plot the frequency response respectly
figure 
plot(wk/pi,20*log10(abs(Hy)),'o');
hold on
plot(W/pi,20*log10(abs(H)),'g');
hold off

%% Problem 6.4 Bandpass Filtering
clear
% load the data
load b3pulses.mat
% calculate the FFT of the input signal and plot it in frequency domain
ffty = fft(x);
plot(1:length(ffty)/2,abs(ffty(1:length(ffty)/2)));
% define the frequency parameters
f1 = [5000  8000];
f2 = [10500 15500];
f3 = [18000 20000];
fs = 80000;
% calculate the bandwidth of each bandpass pulses
delta_w1 = 2*pi*(f1(2)-f1(1))/fs;
delta_w2 = 2*pi*(f2(2)-f2(1))/fs;
delta_w3 = 2*pi*(f3(2)-f3(1))/fs;
% calculate the center frequency of each bandpass pulses
w1 = 2*pi*(f1(2)+f1(1))/(2*fs);
w2 = 2*pi*(f2(2)+f2(1))/(2*fs);
w3 = 2*pi*(f3(2)+f3(1))/(2*fs);
% calculate the bandwidth parameter r for each bandpass pulses
syms r;
r1 = double(solve(2*(1-r)/sqrt(r)-delta_w1));
r2 = double(solve(2*(1-r)/sqrt(r)-delta_w2));
r3 = double(solve(2*(1-r)/sqrt(r)-delta_w3));
% calculate the filter coefficients
b  = [1 0 -1];
a1 = [1 -2*r1*cos(w1) r1^2];
a2 = [1 -2*r2*cos(w2) r2^2];
a3 = [1 -2*r3*cos(w3) r3^2];
% calculae the frequency responses and plot them
[H1,W1] = freqz(b,a1);
[H2,W2] = freqz(b,a2);
[H3,W3] = freqz(b,a3);
figure
plot(W1/pi, 20*log10(abs(H1)),'r--',...
     W2/pi, 20*log10(abs(H2)),'g-.',...
     W3/pi, 20*log10(abs(H3)),'k:'...
    );
% filter the signal x using the filters designed above
y1 = filter(b,a1,x);
y2 = filter(b,a2,x);
y3 = filter(b,a3,x);
ffty1 = fft(y1);
ffty2 = fft(y2);
ffty3 = fft(y3);
% plot the signals in time domain
figure
subplot(2,2,1);
plot(x);
subplot(2,2,2);
plot(y1,'r--');
subplot(2,2,3);
plot(y2,'g-.');
subplot(2,2,4);
plot(y3,'k:');
% plot th signals in frequency domain
figure
subplot(2,2,1);
plot(1:length(ffty)/2,abs(ffty(1:length(ffty)/2)));
subplot(2,2,2);
plot(1:length(ffty1)/2,abs(ffty1(1:length(ffty1)/2)),'r--');
subplot(2,2,3);
plot(1:length(ffty2)/2,abs(ffty2(1:length(ffty2)/2)),'g-.');
subplot(2,2,4);
plot(1:length(ffty3)/2,abs(ffty3(1:length(ffty3)/2)),'k:');




