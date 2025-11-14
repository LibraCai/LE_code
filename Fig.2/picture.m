clc;clear;

n = 2;
N = 100000;
NN = 100;
lambda111 = [];



%parfor ii = 1:NN

%ii
    rho = [];
for i = 1:n
    rho = [rho; 1+i.^(-1)];
end

lambda_init = [];
for i = 1:n
    lambda_init = [lambda_init; 0.5+(1./i)];
end

sigma = 1;
eps = 0.3;
T = 10;

delta_t = T/N;

k = 0.5;

t_s = linspace(0,T, N+1);

lambda = lambda_init;

for i = 1:N
    
    size_current = size(lambda);
    lambda_last = lambda(:,size_current(2));
    lambda_new = lambda_last + F(lambda_last, rho, eps, k).*delta_t + sigma.*((1+2.*eps).^0.5).*sqrt(delta_t).*mvnrnd(zeros(n,1), eye(n), 1)';
    lambda = [lambda lambda_new];
end

%lambda111 = [lambda111; lambda]

%end

lll = max(lambda, [], 1);

figure;
plot(t_s, lll, 'b-', 'DisplayName', 'n = 2','LineWidth',2)
xlabel('t')
ylabel('\lambda_1')
hold on;

%legend('n = 2', 'n = 5', 'n = 7', 'n = 11', 'n = 13');
hold off;

