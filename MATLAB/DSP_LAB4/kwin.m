function [ kw, N, belta ] = kwin( deltaw, A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% A must in db

N = ceil( (A-8) / (2.285*deltaw) + 1 );
if A < 21
   belta = 0; 
elseif A >21 && A < 50
   belta = 0.5842*(A-21)^0.4 + 0.07866*(A-21);
elseif A >= 50
   belta = 0.1102*(A-8.7);
end

n = 0:1:N-1;
kw = besseli(0,belta*sqrt(1-(2*n/(N-1)-1).^2))/besseli(0,belta);
end

