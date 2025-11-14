function Ly_exp = ly_exp_numerical(f_x,sigma,t_max)

%Computing the stochastic Lyapunov exponents for the Langevin system y' = f(y) + σ dW(t).
%σ is the diffusion coefficient, and dW(t) represents Gaussian colored noise.

dt = 0.01;
numTimeSteps = ceil(t_max/dt);


single_path_X_total = [];
single_path_Y_total = []; 
N_trace = 3; % The dimension of the system is 3; therefore, the number of generated trajectories is also 3.
ndim = N_trace;

init_Xs = [0 0 1;
           1 1 1;
           1 0 0];

init_Ys = init_Xs + eye(ndim);

%First, stochastic trajectories for dx(t) = F(x(t)) + σ dW₁(t) are generated,
for n_trace = 1:N_trace

single_path_X = [];
init_X = init_Xs(:,n_trace);

X_p = init_X;

for n_steps = 1:numTimeSteps

dXdt = f_x(X_p);
%dx(t) = F(x(t))dt+σ dW_1(t)
dX = dXdt.*dt + sigma .* (- (1/sqrt(2))*randn(ndim,1)*dt + sqrt(dt)*randn(ndim,1) );

X_p = X_p + dX;

single_path_X = [single_path_X X_p];

end

single_path_X_total = [single_path_X_total; single_path_X];

end


% Next, stochastic trajectories for dy(t) = F(y(t))+σ dW_2(t) are generated,
for n_trace = 1:N_trace

single_path_Y = [];
init_Y = init_Ys(:,n_trace);

Y_p = init_Y;

for n_steps = 1:numTimeSteps

dYdt = f_x(Y_p);
%dy(t) = F(y(t))dt+σ dW_2(t)
dY = dYdt.*dt + sigma .* (- (1/sqrt(2))*randn(ndim,1)*dt + sqrt(dt)*randn(ndim,1) );

Y_p = Y_p + dY;

single_path_Y = [single_path_Y Y_p];

end

single_path_Y_total = [single_path_Y_total; single_path_Y];

end

X_t = single_path_X_total(:,end);
Y_t = single_path_Y_total(:,end);

e_t = Y_t - X_t;

U_t = [e_t(1:3) e_t(4:6) e_t(7:9)];

%H_t = (U_t*U_t')^(1/(2*t_max));
H_t = (U_t*U_t');

eigenvalues = eig(H_t);

Ly_exp = log(eigenvalues)/(2*t_max);