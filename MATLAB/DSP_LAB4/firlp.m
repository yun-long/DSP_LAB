function [ h ] = firlp( N, wc )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% if N is odd
if mod(N,2) ~= 0
    n = -(N-1)/2:(N-1)/2;
% if N is even
else
    n = -N/2:N/2;
end
% generate the truncated impulse response of an ideal FIR filter
w = wc/pi;
h = w*sinc(w*n);
end

