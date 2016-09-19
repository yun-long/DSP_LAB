function [ Cxx ] = spec1( xn, L, win )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = length(xn);
M = N/L;
m = 1:M;
for l = 1:L
    xl(l,m) = xn(m+(l-1)*M);
    winsum  = sum(abs(win).^2);
    winsum  = 1/winsum;
    Xwin(l,m) = abs(fft(win.*xl(l,m))).^2;
    Ixx(l,m) = winsum*Xwin(l,m);
end
Cxx(m) = sum(Ixx(:,m))/L;
end

