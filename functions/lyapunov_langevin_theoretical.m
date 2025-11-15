function Ly_exp = lyapunov_langevin_theoretical(n,f_x,J_x,t_max,ystart,D_diff)

rho = lyapunov_langevin_numerical(n,f_x,J_x,t_max,ystart,0);

%for initial time,all the \mu_i are zero.
dt = 0.005;
numTimeSteps = ceil(t_max/dt);

%single_path_mu = [];
ndim = n;
init_mu = zeros(ndim,1);

mu_p = init_mu;

for n_steps = 1:numTimeSteps
t = n_steps*dt;
dmudt = F_theory(mu_p,rho,D_diff,t);

dmu   = dmudt.*dt + D_diff .* (- (1/sqrt(2))*randn(ndim,1)*dt + sqrt(dt)*randn(ndim,1) );

mu_p = mu_p + dmu;

%single_path_mu = [single_path_mu mu_p];
end

%Ly_exp = single_path_mu(:,end);
Ly_exp = mu_p;
end



