clc;clear;close all;

Tra_grid = 80;
sigma_diff = 0.8;
T_max = 80;
ndim = 3;
y_init = [0.1 0 0.1];
%[~,Res]=lyapunov(ndim,@MyLorenz,@ode45,0,0.01,500,[0 1 0],5);

%Ly_determin=Res(end,:)

%Ly_numerical = ly_exp_numerical(@Lorenz,sigma_diff,T_max);
%Ly_theoretic = ly_exp_theoretical(@Lorenz,sigma_diff,T_max);

%Ly_numerical = sort(Ly_numerical)
%Ly_theoretic = sort(Ly_theoretic)
N_samples = 1000;

Ly_numericals = [];
Ly_theoretics = [];


parfor ii = 1:N_samples

ii
%Ly_numerical = ly_exp_numerical(@Lorenz,sigma_diff,T_max);
Ly_numerical = lyapunov_langevin_numerical(ndim,@Lorenz,@Lorenz_Jacobi,T_max,y_init,sigma_diff);

Ly_numerical = sort(Ly_numerical);

Ly_numericals = [Ly_numericals Ly_numerical];

%Ly_theoretic = ly_exp_theoretical(@Lorenz,sigma_diff,T_max);
Ly_theoretic = lyapunov_langevin_theoretical(ndim,@Lorenz,@Lorenz_Jacobi,T_max,y_init,sigma_diff);

Ly_theoretic = sort(Ly_theoretic);

Ly_theoretics = [Ly_theoretics Ly_theoretic];


end

mu1_min = min([Ly_numericals(1,:) Ly_theoretics(1,:)]);
mu1_max = max([Ly_numericals(1,:) Ly_theoretics(1,:)]);

mu2_min = min([Ly_numericals(2,:) Ly_theoretics(2,:)]);
mu2_max = max([Ly_numericals(2,:) Ly_theoretics(2,:)]);

mu3_min = min([Ly_numericals(3,:) Ly_theoretics(3,:)]);
mu3_max = max([Ly_numericals(3,:) Ly_theoretics(3,:)]);


%%先搞定两种方法对应到的mu1-mu2的联合概率分布

%数值
num_points12_numer = zeros(Tra_grid,Tra_grid);

for ii = 1:N_samples

Ly_numerical = [Ly_numericals(1,ii);
                Ly_numericals(2,ii)];

AA = floor((Ly_numerical(1) - mu1_min)/(mu1_max - mu1_min)*Tra_grid)+1;
BB = floor((Ly_numerical(2) - mu2_min)/(mu2_max - mu2_min)*Tra_grid)+1;

if (AA <= Tra_grid) && (BB <= Tra_grid)
        num_points12_numer(AA,BB) = num_points12_numer(AA,BB) + 1;
end

end

Z_12_numer = sum( num_points12_numer(:) ).*(mu1_max - mu1_min).*(mu2_max - mu2_min)./(Tra_grid^2);
P_12_numer = num_points12_numer ./Z_12_numer;

%理论
num_points12_theo = zeros(Tra_grid,Tra_grid);

for ii = 1:N_samples

Ly_theoretic = [Ly_theoretics(1,ii);
                Ly_theoretics(2,ii)];

AA = floor((Ly_theoretic(1) - mu1_min)/(mu1_max - mu1_min)*Tra_grid)+1;
BB = floor((Ly_theoretic(2) - mu2_min)/(mu2_max - mu2_min)*Tra_grid)+1;

if (AA <= Tra_grid) && (BB <= Tra_grid)
        num_points12_theo(AA,BB) = num_points12_theo(AA,BB) + 1;
end

end

Z_12_theo = sum( num_points12_theo(:) ).*(mu1_max - mu1_min).*(mu2_max - mu2_min)./(Tra_grid^2);
P_12_theo = num_points12_theo ./Z_12_theo;

x1_12 = linspace(mu1_min,mu1_max,Tra_grid+1);
x2_12 = linspace(mu2_min,mu2_max,Tra_grid+1);
[X1_12, X2_12] = meshgrid(x1_12(1:Tra_grid), x2_12(1:Tra_grid));

figure;
subplot(2,2,1);
surf(X1_12, X2_12, P_12_theo,'EdgeColor','none');
xlabel('\mu_1');
ylabel('\mu_2');
zlabel('Prob');
title('\mu_1 -- \mu_2 marginal distribution(theoretical)');
xlim([mu1_min mu1_max]);ylim([mu2_min mu2_max]);

subplot(2,2,2);
surf(X1_12, X2_12, P_12_numer,'EdgeColor','none');
xlabel('\mu_1');
ylabel('\mu_2');
zlabel('Prob');
title('\mu_1 -- \mu_2 marginal distribution(numerical)');
xlim([mu1_min mu1_max]);ylim([mu2_min mu2_max]);

%%下面搞定两种方法对应到的mu2-mu3的联合概率分布

%数值
num_points23_numer = zeros(Tra_grid,Tra_grid);

for ii = 1:N_samples

Ly_numerical = [Ly_numericals(2,ii);
                Ly_numericals(3,ii)];

AA = floor((Ly_numerical(1) - mu2_min)/(mu2_max - mu2_min)*Tra_grid)+1;
BB = floor((Ly_numerical(2) - mu3_min)/(mu3_max - mu3_min)*Tra_grid)+1;

if (AA <= Tra_grid) && (BB <= Tra_grid)
        num_points23_numer(AA,BB) = num_points23_numer(AA,BB) + 1;
end

end

Z_23_numer = sum( num_points23_numer(:) ).*(mu2_max - mu2_min).*(mu3_max - mu3_min)./(Tra_grid^2);
P_23_numer = num_points23_numer ./Z_23_numer;


%理论
num_points23_theo = zeros(Tra_grid,Tra_grid);

for ii = 1:N_samples

Ly_theoretic = [Ly_theoretics(2,ii);
                Ly_theoretics(3,ii)];

AA = floor((Ly_theoretic(1) - mu2_min)/(mu2_max - mu2_min)*Tra_grid)+1;
BB = floor((Ly_theoretic(2) - mu3_min)/(mu3_max - mu3_min)*Tra_grid)+1;

if (AA <= Tra_grid) && (BB <= Tra_grid)
        num_points23_theo(AA,BB) = num_points23_theo(AA,BB) + 1;
end

end

Z_23_theo = sum( num_points23_theo(:) ).*(mu2_max - mu2_min).*(mu3_max - mu3_min)./(Tra_grid^2);
P_23_theo = num_points23_theo ./Z_23_theo;

x1_23 = linspace(mu2_min,mu2_max,Tra_grid+1);
x2_23 = linspace(mu3_min,mu3_max,Tra_grid+1);
[X1_23, X2_23] = meshgrid(x1_23(1:Tra_grid), x2_23(1:Tra_grid));


subplot(2,2,3);
surf(X1_23, X2_23, P_23_theo,'EdgeColor','none');
xlabel('\mu_2');
ylabel('\mu_3');
zlabel('Prob');
title('\mu_2 -- \mu_3 marginal distribution(theoretical)');
xlim([mu2_min mu2_max]);ylim([mu3_min mu3_max]);

subplot(2,2,4);
surf(X1_23, X2_23, P_23_numer,'EdgeColor','none');
xlabel('\mu_2');
ylabel('\mu_3');
zlabel('Prob');
title('\mu_2 -- \mu_3 marginal distribution(numerical)');
xlim([mu2_min mu2_max]);ylim([mu3_min mu3_max]);
