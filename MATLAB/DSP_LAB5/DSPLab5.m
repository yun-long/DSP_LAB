%% DSP Lab 5 IIR Filter Design Using Analog Prototypes!
clear all;
close all;
%% Problem 5.1 Impulse Invariance Method and Bilinear Transformation
% Parameters specification
N  = 4;
wc = 0.3*pi;
Td = 2;
k = 0:N-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step_1: Translate the discrete time filter specifications to the analog
% design domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency transformation
Omegaci = wc/Td;
Omegacb = 2*tan(wc/2)/Td;
% roots of the denominator polynomial
ski = Omegaci*exp(1j*pi*(2*k+N+1)./(2*N));
skb = Omegacb*exp(1j*pi*(2*k+N+1)./(2*N));
% polynomial representation of the denominator
phi = poly(ski);
phb = poly(skb);
% analog filter coefficients
ai = phi;
bi = Omegaci^N;
ab = phb;
bb = Omegacb^N;
% calculate trasfer function in continuous time domain
hti = tf(bi,ai);
htb = tf(bb,ab);
% calculate the frequency response in continuous time domain
[Hti, Wti] = freqs(bi,ai,w);
[Htb, Wtb] = freqs(bb,ab,w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step_2: Select an analog prototype(e.g. Butterworth)
% Step_3: Apply the continuous to discrete-time transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the filter coefficients of discrete time domain 
[bzi, azi] = impinvar(bi,ai,1/Td);
[bzb, azb] = bilinear(bb,ab,1/Td);
bzi = real(bzi);
azi = real(azi);
bzb = real(bzb);
azb = real(azb);
% calculate the transfer function in discrete time domain
% Hzi = tf(bzi,azi,Td);
% Hzb = tf(bzb,azb,Td);
% delta function generation
n = 0:50;
delta0 = +(n==0);
% calculate the impulse response in discrete time domain
hni = filter(bzi,azi,delta0);
hnb = filter(bzb,azb,delta0);
% calculate the frequency respone in discrete time domain
[Hni, Wni] = freqz(bzi,azi);
[Hnb, Wnb] = freqz(bzb,azb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step_4: Display the result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the impulse response h(t) and h(n)
%impulse invariance vs bilinear transform--h(t)
figure
impulse(hti);
hold on
impulse(htb,'r:');
%impulse invariance vs bilinear transform--h(n)
figure
subplot(2,1,1);
stem(n,hni);
subplot(2,1,2);
stem(n,hnb);
% plot the frequency respone |H(jomega)| and |H(jOmega)|
%continuous time domain--|H(jomega)|
figure
subplot(2,1,1);
plot(Wti,20*log10(abs(Hti)),Wtb,20*log10(abs(Htb)),'r:');
subplot(2,1,2);
plot(Wti,angle(Hti),Wti,angle(Htb),'r:');
%discrete time domain--|H(jOmega)|
figure
subplot(2,1,1);
plot(Wni,20*log10(abs(Hni)),Wnb,20*log10(abs(Hnb)),'r:');
subplot(2,1,2);
plot(Wni,angle(Hni),Wnb,angle(Hnb),'r:');

%% Problem 5.2 Characteristics of IIR Filters
clear;
% parameters specification
fsamp = 200;
fp = 32;
ft = 38;
Rp = 1;
Rs = 25;
wp = fp/(fsamp/2);
ws = ft/(fsamp/2);
% determine the order of the butterworth, chebyshev, elliptic filters
[Nbutt, Wnbutt] = buttord(wp,ws,Rp,Rs);
[Ncheb, Wncheb] = cheb1ord(wp,ws,Rp,Rs);
[Nelli, Wnelli] = ellipord(wp,ws,Rp,Rs);
% design the specified filters
[bbutt, abutt] = butter(Nbutt,Wnbutt);
[bcheb, acheb] = cheby1(Ncheb,Rp,Wncheb);
[belli, aelli] = ellip(Nelli,Rp,Rs,Wnelli);
% calculate the frequency response of the filters
[Hwbutt, Wbutt] = freqz(bbutt,abutt);
[Hwcheb, Wcheb] = freqz(bcheb,acheb);
[Hwelli, Welli] = freqz(belli,aelli);
% delta function generation
n = 0:50;
delta0 = +(n==0);
% calculate the impulse response of the filters
hnbutt = filter(bbutt,abutt,delta0);% ???? Some thing wrong with this 
hncheb = filter(bcheb,acheb,delta0);
hnelli = filter(belli,aelli,delta0);
% calculate the transfer function
Hzbutt = tf(bbutt,abutt,1/fsamp);
Hzcheb = tf(bcheb,acheb,1/fsamp);
Hzelli = tf(belli,aelli,1/fsamp);
% display the results
%pole-zero location diagram
figure
pzmap(Hzbutt);
figure
pzmap(Hzcheb);
figure
pzmap(Hzelli);
%impulse response
figure
subplot(3,1,1);
stem(n,hnbutt);
subplot(3,1,2);
stem(n,hncheb);
subplot(3,1,3);
stem(n,hnelli);
%frequency response
figure
subplot(3,1,1);
plot(Wbutt,20*log10(abs(Hwbutt)));
subplot(3,1,2);
plot(Wcheb,20*log10(abs(Hwcheb)));
subplot(3,1,3);
plot(Welli,20*log10(abs(Hwelli)));
% optimum FIR filter desing
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1), 10^(-Rs/20)];
[Npm, fo, ao, Wpm] = firpmord([32, 38],[1,0],dev,fsamp);
bpm = firpm(Npm,fo,ao,Wpm);
figure
freqz(bpm,1,1024,fsamp);
%% Problem 5.3.1 IIR Filtering of Sinusoids
fsamp = 1000;
f1 = 100;
f2 = 150;
n  = 1:300;
wp = 150/(fsamp/2);
ws = 100/(fsamp/2);
rp = 1;
rs = 40;
% determine the order of the filter
[N, wn] = cheb2ord(wp,ws,rp,rs);
% design the filter
[b, a] = cheby2(N,rs,wn,'high');
% calculate the transfer function
Hn = tf(b,a,1/fs);
% time domain
%input signal
xt = 5*sin(2*pi*f1/fsamp*n) + 2*sin(2*pi*f2/fsamp*n);
%filtered signal
yt = filter(b,a,xt);
% frequency domain
Hxt = fft(xt);
N1 = length(Hxt);
w1 = (0:1/N1:1-1/N1)*2*pi;
Hyt = fft(yt);
N2 =length(Hyt);
w2 = (0:1/N2:1-1/N2)*2*pi;
% display the result
%frequency response of the filter
figure
freqz(b,a);
%input signal and output singal in time domain
figure
subplot(2,1,1);
plot(n,xt);
subplot(2,1,2);
plot(n,yt,'r:')
%input signal and output singal in frequency domain
figure
subplot(2,1,1);
plot(w1,abs(Hxt));
subplot(2,1,2);
plot(w2,abs(Hyt),'r:');

%% Problem 5.3.2 IIR Filtering of Sinusoids
N = 10;
fs = 1000;
rp = 1;
rs = 40;
ws = 0.44;
deltaw = 0.08;
w1 = ws - deltaw/2;
w2 = ws + deltaw/2;
% design the filter
[b, a] = ellip(N,rp,rs,[w1, w2],'stop');
% calculate the transfer function
Hz = tf(b,a,1/fs);
% input signal
n = 0:200;
xn = sin(0.44*pi*n);
Hxn = fft(xn);
Nxn = length(Hxn);
wxn = (0:1/Nxn:1-1/Nxn)*2*pi;
% filter the signal
yn = filter(b,a,xn);
Hyn = fft(yn);
Nyn = length(Hyn);
wyn = (0:1/Nyn:1-1/Nyn)*2*pi;
% display the result
%
figure
freqz(b,a);
%
figure
subplot(2,1,1);
plot(n,xn);
subplot(2,1,2);
plot(n,yn,'r:');
%
figure
subplot(2,1,1);
plot(wxn,abs(Hxn));
subplot(2,1,2);
plot(wyn,abs(Hyn),'r:');




