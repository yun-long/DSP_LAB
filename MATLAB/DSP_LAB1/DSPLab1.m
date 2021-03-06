%% Hello this is my first DSP lab. Good luck with it. 

%% Preparation
% % % % % % % % c)
clear
v = -1:0.1:1;
b = sign(v);
for n = 1:length(v)
   if v(n) < 0
       a(n) = -1;
   elseif v(n) > 0
       a(n) = 1;
   elseif v(n) == 0;
       a(n) = 0;
   end
end
% % % % % % % % d)
G = [0.6  1.5  2.3  -0.5;...
     8.2  0.5 -0.1  -2.0;...
     5.7  8.2  9.0   1.5;...
     0.5  0.5  2.4   0.5;...
     1.2 -2.3 -4.5   0.5];
 % of cause it is square
 size(G);
 for n1 = 1:size(G,1)
     for n2 =1:size(G,2)
         if G(n1,n2) == 0.5;
             fprintf('G(%d %d) = 0.5\n', n1, n2);
         end         
         if G(n1,n2) < 0;
             fprintf('G(%d %d) < 0 \n', n1, n2);
         end
     end
 end
 % e)
 help arith;
 help relop;
 help if;
 help for;
 help plot;
%% Problem 4.1 Magic Matrics
clear
M = magic(5);
sum1 = sum(M);
sum2 = sum(M');
row1 = M(1,:);
colum3 = M(:,3); 
columrow = M(2:end,1:3);
i = 1;
j = 1;
m = 1;
n = 1;
for row = 1:size(M,1)
    for col = 1:size(M,2)
        if M(row,col) > 10
           indices1(i,j) = row;
           indices1(i,j+1) = col;
           i = i + 1;
           j = 1;
        end
        if M(row,col) < 4
           indices2(n,m) = row;
           indices2(n,m+1) = col;
           n = n + 1;
           m = 1;
        end
    end
end
%% Problem 4.2 Magic Matrics
clear
N = 12;
fn = fibonacci(N);
for n = 2:N
    ratio(n-1) = fn(n)/fn(n-1); 
end
figure
plot(ratio);
hold on
x = 0:0.1:length(ratio);
Phi = (sqrt(5)+1)/2;
plot(x,Phi,'r.-');
%% Problem 4.3 Statistical Measurements
clear
x     = rand(1000,1);
maxx  = max(x);
minx  = min(x);
meanx = mean(x);
stdx  = std(x);
%
y     = 4*x-2;
maxy  = max(y);
miny  = min(y);
meany = mean(y);
stdy  = std(y);
%
xn    = randn(1000,1);
subplot(2,1,1);
hist(x);
subplot(2,1,2);
hist(xn);
%% Problem 4.4 An Optimization Example
clear
V = 330;
r = 0.5:0.1:10;
A = 2*pi*r.^2+2*V./r;
% find the indices for the minimum of A(r)
n = find(A == min(A));
h = V./(pi*min(A)^2);
plot(r,A);
hold on
plot(r(n),A(n),'r*');
text(r(n),A(n)+20,'minimu');
%% Problem 4.5 The Moving Average
clear
load glob_warm.mat
m = 7;
Xn = movingAverage(Ta,m);
figure
plot(Ta,'b--');
hold on
plot(Xn,'r.-');
hold off
%% Problem 4.6 Signal Processing Example
clear
n = 0:100;
F = 1;
T = 0.05;
% s = sin(2*pi*f*t) = sin(2*pi*f*n*T)
s = sin(2*pi*F*n*T); 
% time domain
figure
subplot(3,1,1);
plot(n,s);
subplot(3,1,2);
plot(n*T,s);
subplot(3,1,3);
stem(n,s);
% frequency domain
S = fft(s,128);
P = S.*conj(S);
w = (0:127)/128;
figure
subplot(2,1,1);
plot(2*w,P);
subplot(2,1,2);
plot(w/T,P);
% adding noise
s2 = s + sin(2*pi*4*n*T);
figure
plot(n,s);
hold on
plot(n,s2,'r.-');
hold off
% filter design
b = [1 1 1 1]/4;
a = 1;
[H,w1] = freqz(b,a);
figure
plot(w1/(2*pi*T),abs(H));
% apply filter to signal
sf = filter(b,a,s2);
figure 
plot(n,sf);
hold on
plot(n,s,'r.-');
hold off














