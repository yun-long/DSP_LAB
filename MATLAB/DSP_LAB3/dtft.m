function [ y, w ] = dtft( x, m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
w = 2*pi*(0:1/m:1-1/m);
y = fft(x,m);
end

