% Set parameters for forward/inverse problems
params_fwd = set_parameters(1,2,1,1,1,0.42,0.25,0.1,1.5,1000,1000,20,[0.02, 0.05, 0.125, 0.245, 0.405]);
params_inv = set_parameters(1,2,1,1,1,0.42,0.25,0.1,1.5,250,250,20,[0.02, 0.05, 0.125, 0.245, 0.405]);

% Set true permeability to be recovered
u_true = sim_permeability(zeros(1,params_fwd.Nx),set_covariance_matrix(params_fwd)); 
plot(params_fwd.x_locations,u_true)
xlim([0,1])
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
disp(ups_true(end))

% Perform LMAP
[u_map,C_map,timer_LMAP,iterations_LMAP] = run_LMAP(u_bar,C,params_inv,experiment,1);

% Perform EKI
[U_EKI,timer_EKI,iterations_EKI] = run_EKI(u_bar,C,params_inv,experiment,1);

% Perform RML
[U_RML,timer_RML,iterations_RML] = run_RML(u_bar,C,params_inv,experiment,1);

% Perform MCMC
[U_MCMC,timer_MCMC,iterations_MCMC] = run_MCMC(u_bar,C,params_inv,experiment,1);

% Plot together
comparison_plotter(u_map,C_map,U_EKI,U_RML,U_MCMC,params_inv,experiment)
exportgraphics(gcf,'C:\Users\pmymc12\OneDrive - The University of Nottingham\PhD\LMAP Paper\Figures\comparison1D.eps')


% Perform sequential LMAP
[u_LMAPs,std_LMAPs,iterations_LMAPs,timer_LMAPs] = run_LMAP_sequential(u_bar,C,params_inv,experiment,1);
exportgraphics(gcf,'C:\Users\pmymc12\OneDrive - The University of Nottingham\PhD\LMAP Paper\Figures\sequentialLMAP1D.eps')