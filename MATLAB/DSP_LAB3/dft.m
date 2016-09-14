function [ y ] = dft( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = length(x);
W = exp(-1j*2*pi/N);
X = zeros(N,N);
for n1 = 1:N
   for n2 = 1:N
       X(n1,n2) = W^((n1-1)*(n2-1));
   end
end
y = X*x';
end

