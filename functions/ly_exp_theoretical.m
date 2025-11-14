function Ly_exp = ly_exp_theoretical(f_x,sigma,t_max)

rho = ly_exp_numerical(f_x,0,t_max);

%for initial time,all the \mu_i are zero.
dt = 0.01;
numTimeSteps = ceil(t_max/dt);

single_path_mu = [];
ndim = 3;
init_mu = zeros(ndim,1);

mu_p = init_mu;

for n_steps = 1:numTimeSteps
t = n_steps*dt;
dmudt = F_theory(mu_p,rho,sigma,t);

dmu   = dmudt.*dt + sigma .* (- (1/sqrt(2))*randn(ndim,1)*dt + sqrt(dt)*randn(ndim,1) );

mu_p = mu_p + dmu;

single_path_mu = [single_path_mu mu_p];
end

Ly_exp = single_path_mu(:,end);
end



