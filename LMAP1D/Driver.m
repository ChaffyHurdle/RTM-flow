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

% Set true permeability to be recovered
u_true = sim_permeability(zeros(1,params_fwd.Nx)',set_covariance_matrix(params_fwd));
figure(1)
plot(params_fwd.x_locations,u_true)
xlim([0,params_fwd.L])
ylim([-1.5,1.5])

% Define mean function and covariance matrix
C = set_covariance_matrix(params_inv);
u_bar = zeros(1,params_inv.Nx);
params_inv.C = C;
params_inv.u_bar = u_bar;

% Generate experimental data and save in object
[p_true,ups_true] = forward_map(u_true,params_fwd);
[d,Sigma] = set_data(p_true,0.00,0.01,params_fwd);
experiment.d = d;
experiment.Sigma = Sigma;
experiment.u_true = u_true;
experiment.params_fwd = params_fwd;
experiment.ups_true = ups_true;

%%
% Perform sequential LMAP
[u_LMAPs,std_LMAPs,iterations_LMAPs,timer_LMAPs] = run_LMAP_sequential(u_bar,C,params_inv,experiment,1);
exportgraphics(gcf,'C:\Users\pmymc12\OneDrive - The University of Nottingham\PhD\LMAP Paper\Figures\sequentialLMAP1D.eps')

%% Perform LMAP
[u_map,C_map,timer_LMAP,iterations_LMAP] = run_LMAP(u_bar,C,params_inv,experiment,0);

%% Perform EKI
[U_EKI500,timer_EKI500,iterations_EKI500] = run_EKI(u_bar,C,params_inv,experiment,500,1);
[U_EKI1000,timer_EKI1000,iterations_EKI1000] = run_EKI(u_bar,C,params_inv,experiment,1000,1);
[U_EKI5000,timer_EKI5000,iterations_EKI5000] = run_EKI(u_bar,C,params_inv,experiment,5000,1);

%% Perform RML
[U_RML500,timer_RML500,iterations_RML500] = run_RML(u_bar,C,params_inv,experiment,500,1);
[U_RML1000,timer_RML1000,iterations_RML1000] = run_RML(u_bar,C,params_inv,experiment,1000,1);
[U_RML5000,timer_RML5000,iterations_RML5000] = run_RML(u_bar,C,params_inv,experiment,5000,1);

%% Perform MCMC
[U_MCMC,timer_MCMC,iterations_MCMC] = run_MCMC(u_bar,C,params_inv,experiment,1);

%% Plot together
comparison_plotter(u_map,C_map,U_EKI5000,U_RML1000,U_MCMC,params_inv,experiment)
exportgraphics(gcf,'C:\Users\pmymc12\OneDrive - The University of Nottingham\PhD\LMAP Paper\Figures\comparison1D.eps')

disp('Means')
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI500,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI1000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_EKI5000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_RML500,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_RML1000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )
disp( sqrt(params_inv.dx*sum( (mean(U_MCMC,2) - mean(U_RML5000,2)).^2 ))/sqrt(params_inv.dx*sum( (mean(U_MCMC,2)).^2 )) )

disp('Stds')
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI500,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI1000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_EKI5000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_RML500,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_RML1000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
disp( sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2)) - sqrt(var(U_RML5000,0,2))).^2 ))/sqrt(params_inv.dx*sum( (sqrt(var(U_MCMC,0,2))).^2 )) )
