function F_val = F_theory(mu,rho,sigma,t)
epss = 10.^(-6);
F_val = -mu+rho;

[ndim,~] = size(mu);

for i = 1:ndim

for j = 1:ndim

if ~(j == i)
F_val(i) = F_val(i) + ...
    (1/(2*sqrt(2)))*(sigma.^2) ...
    /tanh(t.* mu(i) - t.* mu(j) + epss);

end

end
end

end