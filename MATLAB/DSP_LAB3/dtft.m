function [ y, w ] = dtft( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = length(x);
m = 16*N;
w = 2*pi*(0:1/m:1-1/m);
y = fft(x,m);
end

