function [ w ] = convmat(u,v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


u0 = zeros(1,length(v)-1);
u1 = [u u0];
for n = 1:length(v)
   u2(:,n) = circshift(u1',n-1); 
end
w  = u2*v';
end

