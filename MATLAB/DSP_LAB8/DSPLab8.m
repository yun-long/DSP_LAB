%% LAB8 Parametric Spectrum Estimation

%% Problem 6.1 AR Process
N = 1024;
Z = randn(N,1);
b = 1;
a = [1 0.5 0.7 0.2];
X = filter(b,a,Z);
p = 3;
% estimate the AR parameters using build in AR-Yule-Walker function
aco = aryule(X,3);
% 
[Pxx, F] = pyulear(X,3);
figure
plot(F/(pi),10*log10(Pxx));
%true spectrum
[H1, w1] = freqz(b,a);
figure
plot(w1/pi,20*log10(abs(H1)));
%estimated spectrum
[H2, w2] = freqz(b,aco);
figure
plot(w2/pi,20*log10(abs(H2)));


%% Problem 6.2 MA Process
load ma3.mat;
%parameters specification
a = 1;
b = [1, 0.4, -0.2, 0.15];
q = 3;
% 
[H, F] = freqz(b,a);
figure
plot(F,20*log10(abs(H)));
hold on
% estimate the parameter
be = maparam(x,q);
[H1,F1] = freqz(be,a);
plot(F1,20*log10(abs(H1)),'r:');

%% Problem 6.3 Order Selection
load arunknown.mat
M = 10;
N = length(X);

for m = 1:M
   ak = aryule(X,m);   
   sum(m) = 0;
   for n = m:N-1
      sumak = 0;
      for k = 1:m
         sumak = sumak + ak(k)*X(n+1-k); 
      end
      sum(m) = sum(m) + (X(n)+sumak).^2;
   end
   deltamsqrt(m) = sum(m)/(N-m);
   AIC(m) = log(deltamsqrt(m)) + m*2/N;
   MDL(m) = log(deltamsqrt(m)) + m*log(N)/N;
end
m = 1:M;
figure
plot(m,AIC,m,MDL,'r:');


%% Problem6.4 Speech Signal Analysis
clc
clear
load s5.mat
SH = s5(15600:16300);
AA = s5(16800:17500);
%
p = 10;
b = 1;
aSH = aryule(SH,p);
aAA = aryule(AA,p);
[H1, w1] = freqz(b,aSH);
[H2, w2] = freqz(b,aAA);
[pxxSH, wSH] = pwelch(SH);
[pxxAA, wAA] = pwelch(AA);
% estimate the order
M = 10;
[AICSH, MDLSH] = AIC_MDL(SH,M);
[AICAA, MDLAA] = AIC_MDL(AA,M);
[min1, pSH] = min(AICSH);
[min2, pAA] = min(AICAA);
aSH2 = aryule(SH,pSH);
aAA2 = aryule(AA,pAA);
[H3, w3] = freqz(b,aSH2);
[H4, w4] = freqz(b,aAA2);
figure
plot(w3,20*log10(abs(H3)),w1,20*log10(abs(H1)),'black--',wSH,10*log10(pxxSH),'r:');
figure
plot(w4,20*log10(abs(H4)),w2,20*log10(abs(H2)),'black--',wAA,10*log10(pxxAA),'r:');

