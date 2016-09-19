function [ Cxx ] = spec2( xn, m)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N   = length(xn);
Ixx = fft(xn);
Ixx = abs(Ixx).^2/N;

for n = 1:N
    if (((n-m)<1) || ((n+m)>N) )
        Cxx(n) = Ixx(n);
    else
        Cxx(n) = sum(Ixx(n-m:n+m))/(2*m+1);
    end
    
end

end

