function [ Fn ] = fibonacci( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fn(1) = 1;
Fn(2) = 1;
for n = 3:N;
    Fn(n) = Fn(n-1) + Fn(n-2);
end
end

