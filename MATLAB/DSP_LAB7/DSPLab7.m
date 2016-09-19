%% DSP LAB 7 Non-Parametric Spectrum Estimation


%% Problem 4.1 Windowed, Averaged Periodograms
load ar7.mat;
N = length(X);
L = 5;
M = N/L;
% rectangular window and Hamming Windows
Rwin = rectwin(M)';
Hwin = hamming(M)';
CxxRwin = spec1(X,L,Rwin);
CxxHwin = spec1(X,L,Hwin);
figure
subplot(2,1,1);
plot(CxxRwin);
subplot(2,1,2);
plot(CxxHwin);
% 
n = norminv([0.025 0.975],0,1);
log1 = 10*log10(CxxHwin)-10*log(n(1)/sqrt(L));
log2 = 10*log10(CxxHwin)+10*log(n(2)/sqrt(L));
w = (0:1/M:1-1/M)*2*pi;
figure
plot(w,10*log10(CxxHwin),w,log1,'r:',w,log2,'black');

%% Problem 4.2 Smoothed Periodograms
load ar7.mat;
m = 20;
Cxx = spec2(X,m);
figure
plot(Cxx);

%% Probelm 4.3 Sinusoids in Noise 






    





