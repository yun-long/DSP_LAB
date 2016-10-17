function [ GCC ] = genxcorr( x1,x2,lag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = length(x1);
% % number of segments
% M = 8;
% m0 = 1:M;
% % length of segments
% L = N/M;
% l0 = 1:L;
% % fourier transform of signal x1 and x2;
% for m = 1:M
%    xm1(m,:) = x1(l0+(m-1)*L);
%    xm2(m,:) = x2(l0+(m-1)*L);
%    FFTx1 = fft(xm1(m,:),L);
%    FFTx2 = fft(xm2(m,:),L);
%    Ix1x2(m,:) = (FFTx1.*conj(FFTx2))/L;
% end
epsilon = 0.0000007;

FFTx1 = fft(x1,N);
FFTx2 = fft(x2,N);
Ix1x2 = (FFTx1.*conj(FFTx2))/N;

gxx = ifft(Ix1x2./(abs(Ix1x2)+epsilon));

gxx = ifftshift(gxx);

figure
plot(gxx);

GCC = gxx(N/2-lag:N/2+lag);
end

