clc;clear;close all;
NNN = 30;
T_maxs = linspace(30,300,NNN);
N_tmaxs = length(T_maxs);
mean_theoretics = [];
mean_numericals = [];
y_init = [1 0 1];
ndim = 3;
for i_tmax = 1:N_tmaxs

%Tra_grid = 80;
sigma_diff = 0.8;
%sigma_diff = 0.005;
T_max = T_maxs(i_tmax);
ndim = 3;
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



mean_theoretics = [mean_theoretics mean(Ly_theoretics,2)];

mean_numericals = [mean_numericals mean(Ly_numericals,2)];

end

figure;
plot(T_maxs,mean_theoretics(1,:),'r');

hold on;

plot(T_maxs,mean_theoretics(2,:),'g');

hold on;

plot(T_maxs,mean_theoretics(3,:),'b');

hold on;

scatter(T_maxs,mean_numericals(1,:),"square");

hold on;

scatter(T_maxs,mean_numericals(2,:),"diamond");

hold on;

scatter(T_maxs,mean_numericals(3,:),"^");
legend('\mu_1 (theoretical)','\mu_2 (theoretical)','\mu_3 (theoretical)','\mu_1 (numerical)','\mu_2 (numerical)','\mu_3 (numerical)');
xlabel('t');ylabel('\mu');
hold off;
