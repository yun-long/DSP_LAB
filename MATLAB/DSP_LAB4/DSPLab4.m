%% Lab 4: Digital FIR Filter Design

%% Problem 7.1 Windowing Method
% parameters specification
wc = 0.5*pi;
N  = [17, 51];
% generate the truncated impulse response of FIR filter, for different N
h1 = firlp(N(1),wc);
h2 = firlp(N(2),wc);
% generate windows
rectWin1 = rectwin(N(1));
rectWin2 = rectwin(N(2));
bartWin1 = bartlett(N(1));
bartWin2 = bartlett(N(2));
hammWin1 = hamming(N(1));
hammWin2 = hamming(N(2));
% multiply the FIRLp with different windows
hdrect1 = fft(h1.*rectWin1',512);
hdrect2 = fft(h2.*rectWin2',512);
hdbart1 = fft(h1.*bartWin1',512);
hdbart2 = fft(h2.*bartWin2',512);
hdhamm1 = fft(h1.*hammWin1',512);
hdhamm2 = fft(h2.*hammWin2',512);
w = (0:1/512:1-1/512)*2*pi;
% display the frequency responses of FIR low pass filter 
y1 = fft(h1,512);
y2 = fft(h2,512);
w1 = (0:1/length(y1):1-1/length(y1))*2*pi;
w2 = (0:1/length(y2):1-1/length(y2))*2*pi;
figure
plot(w1,abs(y1),w2,abs(y2),'r--');
% find the 3db point
yhdrect1 = 20*log10(abs(hdhamm1));
m = 3;
[Min1, index1] = min(abs(yhdrect1)-m);
% compare the frequency response of different windows with same length
figure
plot(w,20*log10(abs(hdrect1)),'r--',...
     w,20*log10(abs(hdbart1)),'g:',...
     w,20*log10(abs(hdhamm1)),...
     w(index1),yhdrect1(index1),'ro');
% compare the frequency response of the same window in different length
%rectangular window
figure
plot(w,20*log10(abs(hdrect1)),'r--',...
     w,20*log10(abs(hdrect2)));
%bartlett window
figure
plot(w,20*log10(abs(hdbart1)),'r--',...
     w,20*log10(abs(hdbart2)));
%hamming window
figure
plot(w,20*log10(abs(hdhamm1)),'r--',...
     w,20*log10(abs(hdhamm2)));

%% Problem 7.2 Kaiser Window
% parameters definition
%low pass parameters
wp = 0.12*pi;
ws = 0.18*pi;
A  = 30;
deltawlp = ws - wp;
wc = (ws+wp)/2;
%high pass parameters
wch = 0.15*pi;
%bandpass parameters
w0 = 0.25*pi;
deltawbp = 0.3*pi;
wl = w0 - deltawbp/2;
wh = w0 + deltawbp/2;
% create kaiser window
[kaiserWin, N] = kwin(deltawlp, A);
% time scale
n = 0:N-1;
% design a low pass filter
lp = firlp(N,wc);
kaiserlp = lp.*kaiserWin;
% design an all pass filter
if mod(N,2)==0
    del = +(n == N/2);
else
    del = +(n == (N-1)/2);
end
% design a high pass filter
hp = del - lp;
kaiserhp = hp.*kaiserWin;
% design a bandpass filter;
lp1 = firlp(N,wl);
lp2 = firlp(N,wh);
bp  = lp2 - lp1;
kaiserbp = bp.*kaiserWin;
% for displayment
wvtool(kaiserWin);
wvtool(kaiserlp);
wvtool(kaiserhp);
wvtool(kaiserbp);
%% Problem 7.3 Frequency Sampling Method
wc = 0.15*pi;
N = 53;
alpha = (N-1)/2;
x = ceil((wc/pi)*alpha);
%
Hmag = [ones(1,x), zeros(1,alpha+1-x)];
wk = 2*pi/N*(0:alpha);
Hphi = exp(-1i*wk*alpha);
H = Hmag.*Hphi;
H = [H, conj(H(end:-1:2))];
h = real(ifft(H));
%
wvtool(h);
%% Problem 7.4 Parks-McClellan Algorithm
% parameters specification
f = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 0.9 1.0];
a = [0 0   1   1   0   0   0.5 0.5 0   0  ];
N = [10, 25, 50];
% filter design for different filter oder N
h1 = firpm(N(1), f, a);
h2 = firpm(N(2), f, a);
h3 = firpm(N(3), f, a);
% plot the designed filter
wvtool(h1,h2,h3);
%% Problem 7.5 Hands-on Example: Outdoor Recording
load noise_5insx.mat
load speech.mat
x_record = x_speech + 0.5*x_5insx;
% soundsc(x_record);
h = fft(x_record);
N =length(h);
w = (0:1/N:1-1/N)*2;
plot(w,abs(h));
% low pass filter design use window method
wc = 0.5*pi;
hammWin = hamming(N);
hx = firlp(N,wc);
hd = hx.*hammWin';
% filter signal
yout = conv(hd,x_record);
figure
plot(abs(fft(yout)));
soundsc(yout);