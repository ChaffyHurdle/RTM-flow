%% Initialise

% Set parameters for forward/inverse problems
L = 1;
p_I = 2;
p_0 = 1;
mu = 1;
phi = 1;
T = L^2*mu*phi/(2*(p_I-p_0)); % based on prior mean
var_matern = 0.25;
l = 0.1;
nu = 1.5;
nsensors = 20;
nobtimes = 5;
ob_times = ((1:nobtimes)/nobtimes*L).^2*mu*phi/(2*(p_I-p_0));

params_fwd = set_parameters(L,p_I,p_0,mu,phi,T,var_matern,l,nu, ...
    1000,1000,nsensors,ob_times);
params_inv = set_parameters(L,p_I,p_0,mu,phi,T,var_matern,l,nu, ...
    250,250,nsensors,ob_times);

% Sample or load true permeability to be recovered
% u_true = sim_permeability(zeros(1,params_fwd.Nx)',set_covariance_matrix(params_fwd));
u_true = load("data1D/u_true.mat").u_true;
figure(1)
plot(params_fwd.x_locations,u_true)
xlim([0,params_fwd.L])
ylim([-1.5,1.5])

% Define mean function and covariance matrix
C = set_covariance_matrix(params_inv);
u_bar = zeros(1,params_inv.Nx);
params_inv.C = C;
params_inv.u_bar = u_bar;

% Generate experimental data
[p_true,ups_true] = forward_map(u_true,params_fwd);
[d,Sigma] = set_data(p_true,0.00,0.01,params_fwd);
Experiment.d = d;
Experiment.Sigma = Sigma;
Experiment.u_true = u_true;
Experiment.params_fwd = params_fwd;
Experiment.ups_true = ups_true;

% or load previous data
% Experiment = load("data1D/Experiment.mat");

%% Sequential LMAP plot

% Set Levenberg-Marquardt parameters
LevenbergMarquardt.tol_J = 0.01;
LevenbergMarquardt.tol_U = 0.01;
LevenbergMarquardt.N_iter = 100;
LevenbergMarquardt.scaling = 2;
LevenbergMarquardt.alpha0 = 1e4;

% Perform sequential LMAP
[u_LMAPs,std_LMAPs,iterations_LMAPs,timer_LMAPs] = run_LMAP_sequential(u_bar, C, params_inv, Experiment, LevenbergMarquardt, 1);
exportgraphics(gcf,'plots1D\sequentialLMAP1D.eps')
exportgraphics(gcf,'plots1D\sequentialLMAP1D.png')
save('data1D/LMAP_seq_data','u_LMAPs','std_LMAPs','iterations_LMAPs','timer_LMAPs')

%% Perform LMAP
[u_map,C_map,timer_LMAP,iterations_LMAP] = run_LMAP(u_bar,C,params_inv,Experiment,LevenbergMarquardt,0);

%% Perform EKI
[U_EKI500,timer_EKI500,iterations_EKI500] = run_EKI(u_bar,C,params_inv,Experiment,500,1);
[U_EKI1000,timer_EKI1000,iterations_EKI1000] = run_EKI(u_bar,C,params_inv,Experiment,1000,1);
[U_EKI5000,timer_EKI5000,iterations_EKI5000] = run_EKI(u_bar,C,params_inv,Experiment,5000,1);
save('data1D/EKI_data', 'U_EKI500','timer_EKI500','iterations_EKI500', ...
                        'U_EKI1000','timer_EKI1000','iterations_EKI1000', ...
                        'U_EKI5000','timer_EKI5000','iterations_EKI5000')

%% Perform RML
[U_RML1000,timer_RML1000,iterations_RML1000] = run_RML(u_bar,C,params_inv,Experiment,LevenbergMarquardt,1000,1);
save('data1D/RML_data', 'U_RML1000', 'timer_RML1000', 'iterations_RML1000')

%% Perform MCMC
[U_MCMC,timer_MCMC,iterations_MCMC] = run_MCMC(u_bar,C,params_inv,Experiment,1);
save('data1D/MCMC_data', 'U_MCMC', 'timer_MCMC', 'iterations_MCMC')

%% Comparison plot
comparison_plotter(u_map,C_map,U_EKI5000,U_RML1000,U_MCMC,params_inv,Experiment)
exportgraphics(gcf,'plots1D/comparison1D.eps')
exportgraphics(gcf,'plots1D/comparison1D.png')

%% Compute relative errors
rel_error_mean.LMAP = sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - u_map).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 ));
rel_error_mean.EKI500 = sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI500,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 ));
rel_error_mean.EKI1000 = sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI1000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 ));
rel_error_mean.EKI5000 = sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI5000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 ));
rel_error_mean.RML = sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_RML1000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 ));

rel_error_std.LMAP = sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(diag(C_map))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 ));
rel_error_std.EKI500 = sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI500,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 ));
rel_error_std.EKI1000 = sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI1000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 ));
rel_error_std.EKI5000 = sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI5000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 ));
rel_error_std.RML = sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_RML1000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 ));