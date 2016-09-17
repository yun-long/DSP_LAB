function f = firsample( samples )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N = 2*length(samples)-1;
H_d = zeros(1,N);

Phi = pi*(N-1)/N;
H_d(1) = samples(1);

for j = 2:N/2-1
Phase = exp(-1i*(j-1)*Phi);
H_d(j) = samples(j)*Phase;
H_d(N+2-j) = samples(j)*conj(Phase);
end
plot(abs(H_d));
f = real(ifft(H_d));
end

