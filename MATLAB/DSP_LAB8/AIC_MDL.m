function [ AIC, MDL] = AIC_MDL( X,M )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

N = length(X);

for m = 1:M
   ak = aryule(X,m);   
   sum(m) = 0;
   for n = m:N-1
      sumak = 0;
      for k = 1:m
         sumak = sumak + ak(k)*X(n+1-k); 
      end
      sum(m) = sum(m) + (X(n)+sumak).^2;
   end
   deltamsqrt(m) = sum(m)/(N-m);
   AIC(m) = log(deltamsqrt(m)) + m*2/N;
   MDL(m) = log(deltamsqrt(m)) + m*log(N)/N;
end

end

