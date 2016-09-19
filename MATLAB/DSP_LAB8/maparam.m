function [ b ] = maparam( x,q )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Cxx = xcorr(x,q);
rxx = roots(Cxx);
N = length(rxx);
m = 0;
for n = 1:N
   if (abs(rxx(n)) <= 1)
       zz(m+1) = rxx(n);
   end
end
b = 1;
end

