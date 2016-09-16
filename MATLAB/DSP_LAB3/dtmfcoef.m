function [ key ] = dtmfcoef( x, fs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

y = dft(x);
N = length(y);
yAbs = abs(y);
keys = ...
    ['1','2','3','A';
     '4','5','6','B';
     '7','8','9','C';
     '*','0','#','D'];
 
cols = [1209, 1336, 1477, 1633];
rows = [ 697,  770,  852,  941];

f = (0:1/N:1-1/N)*fs;
index = find(yAbs > max(yAbs)/2);
fr = f(index(1));
fc = f(index(2));
for n = 1:4
    if (fr > (rows(n)-60) && fr < (rows(n)+60))
        r = n;
    end
    
    if (fc > (cols(n)-36) && fc < (cols(n)+ 36))
        c = n;
    end
end
key = keys(r,c);
end

