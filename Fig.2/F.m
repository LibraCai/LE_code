function F_val = F(lambda, rho, eps, k)
epsilon = 10^(-10);
n = length(lambda);

F0 = [];
for i = 1:n
    F0 = [F0 lambda];
end
F0 = F0 - lambda';

for i = 1:n
    F0(i,i) = 1;
end
F0 = 1./tanh(F0+epsilon);
for i = 1:n
    F0(i,i) = 0;
end

F_val = -rho + 0.5.*(1+eps).*(k.^2).*sum(F0, 2);
end