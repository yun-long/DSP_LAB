%% DSP LAB10 Localization of Acoustical Sources

%% Problem 3.1 Time-Delay Estimation using the Cross-SOMF
clear
load group04_static.mat
XSpeech = x;
load rec_white.mat
Xnoise = x;
f = 240000/5;
% soundsc(x(:,2),240000/5);
%
% figure
% spectrogram(x(:,1));
index = find(XSpeech(:,1) == max(XSpeech(:,1)));
indexrange = [index-512*(1/2):index+512*(1/2)-1];
s1 = XSpeech(indexrange,1);
s2 = XSpeech(indexrange,2);
%
[Sx1x2,lags]= xcorr(s1,s2);
Sx1x2 = Sx1x2(lags>=0);
lags = lags(lags>=0);
[max, delN] = max(Sx1x2);
delt = 1/f*delN;
figure
plot(lags/512,Sx1x2);
figure
plot(XSpeech(:,1));
figure
plot(s1);
hold on
plot(s2,'r:')

n1 = Xnoise(:,1);
n2 = Xnoise(:,2);
[Sn1n2, lagsn] = xcorr(n1,n2);
figure
plot(Sn1n2);
Sn1n2 = Sn1n2(lags>=-10);
figure
plot(Sn1n2);

%% Problem 3.2 Improvements by GCC
% speech signal: s1, s2
% noise  singal: n1, n2
clear
load group04_static.mat
XSpeech = x;
index = find(XSpeech(:,1) == max(XSpeech(:,1)));
indexrange = [index-512*(0/2):index+512*(2/2)-1];
s1 = XSpeech(indexrange,1);
s2 = XSpeech(indexrange,4);

[Sx1x2,lags]= xcorr(s1,s2);
Sx1x2 = Sx1x2(lags>=0);
lags = lags(lags>=0);
[max, delN] = max(Sx1x2);

GCC = genxcorr(s1,s2,delN);

% figure
% plot(GCC,'r');





