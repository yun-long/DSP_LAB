function [ Xn ] = movingAverage(x, m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(x);
%  
for n = 1:1:N
    if n < m + 1;
        Xn(n) = x(n);
    elseif n > N-m
        Xn(n) = x(n);
    else
        Xn(n) = sum(x([n-m:n+m]))/(2*m+1);
    end
end

end
