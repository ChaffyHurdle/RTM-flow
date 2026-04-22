%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
parpool('Threads');
addpath('Utils')

%% Set case study number and load meshes and sensor locs
case_study = 1; % 1: square, 2: fork, 3: annulus

% Add all relevant meshes and sensor locations as variables
folder_path = strcat('Case',num2str(case_study));
mesh_sensor_path = strcat(folder_path,'/meshes_sensors');
files = dir(fullfile(mesh_sensor_path, '*.mat'));
for k = 1:length(files)
    file_path = fullfile(mesh_sensor_path, files(k).name);
    [~, name, ~] = fileparts(files(k).name);
    tmp = load(file_path);
    var_names = fieldnames(tmp);
    assignin('base', name, tmp.(var_names{1}));
end

% Load inlet/outlet functions
[inlet_func,vent_func] = load_inlets_vents(case_study);

%% Mesh set up (avoiding inverse crimes)
my_forward_mesh = DelaunayMesh(p_fwd,e_fwd,t_fwd);
my_inverse_mesh = DelaunayMesh(p_inv,e_inv,t_inv);

%% Load parameters (prior and Darcy flow). These can also be set manually.
[matern_args,mu,phi,thickness,p_I,p_0,observation_times,T] = load_parameters(case_study);

%% Inverse problem set up (~10 secs)
my_inverse = Inversion(my_forward_mesh,my_inverse_mesh,matern_args);

%% Generate new true permeability, or preload example

% Preload example
K_mat = load(strcat(folder_path,'/Example1/K_true.mat'));
my_inverse.u_true = log(K_mat.K_true);
my_inverse.plot_u_true();

% % Generate new permeability
% while true
%     my_inverse = my_inverse.generate_u();
%     my_inverse.plot_u_true();
%     key_input = input("Confirm ('y') or sample again (any other key)");
%     if key_input == 'y'
%         break
%     end
% end

%% Instantiate physics and pressure classes
sensor_locs = sensor_locs_25; % set to personal preference
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, ...
    p_I, p_0, inlet_func, vent_func, ...
    K_true, sensor_locs, observation_times,T);
my_pressure = Pressure(my_forward_mesh,my_darcy);

%% Plot experimental setup
figure(1)
my_inverse.plot_u_true()
hold on
scatter(sensor_locs(:,1),sensor_locs(:,2),'white','filled')
xlim([min(my_forward_mesh.nodes(:,1)),max(my_forward_mesh.nodes(:,1))])
ylim([min(my_forward_mesh.nodes(:,2)),max(my_forward_mesh.nodes(:,2))])
hold off

%% RTM true simulation (on fine forward mesh)
adjoint = 1; % Set to 1 if using LMAP, 0 if using EKI/MCMC
true_RTMflow = RTMFlow(my_forward_mesh,my_darcy,my_pressure,adjoint);
true_RTMflow = true_RTMflow.run(inf);

%% Generate random data for inverse problem
noise_level = 0.005;
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,noise_level);

%% Perform LMAP (all times)
alpha0 = 1e3; scale = 5; tol_U = 0.025; tol_J = 0.025;
my_lmap = LMAP(my_inverse,my_darcy,alpha0,scale,tol_U,tol_J);
my_lmap = my_lmap.run();

plot_LMAP_seq(my_forward_mesh, my_inverse_mesh, my_darcy, true_RTMflow, my_lmap)

%% Save data
lmap_means = my_lmap.umap_seq;
lmap_vars = zeros(size(lmap_means));
lmap_times = my_lmap.timer_seq;
for j = 1:5
    lmap_vars(:,j) = diag(my_lmap.Cmap_seq(:,:,j));
end
for j = 1:5
    t_index = find(true_RTMflow.times > my_darcy.observation_times(j),1)-1;
    edge_data_t = true_RTMflow.edge_data{t_index};
    nodes_t = edge_data_t(:,2:3);
    save(strcat(folder_path,'/Example2/front',num2str(j),'.mat'),"nodes_t")
end
save(strcat(folder_path,"/Example2/K_true.mat"),'K_true')
save(strcat(folder_path,"/Example2/lmap_means.mat"),'lmap_means')
save(strcat(folder_path,"/Example2/lmap_vars.mat"),'lmap_vars')
save(strcat(folder_path,"/Example2/lmap_times.mat"),'lmap_times')

%% Plot push forward example
u_samples = mvnrnd(my_lmap.umap_seq(:,end),my_lmap.Cmap_seq(:,:,end),1000);
[pressures,flow_fronts] = push_forward(u_samples, my_darcy, my_inverse_mesh);
perturbed_pressures = pressures + normrnd(0,1,size(pressures)).*sqrt(my_inverse.Sigma(:)');
plot_push_forward(perturbed_pressures, flow_fronts, my_darcy, my_inverse_mesh, true_RTMflow)

%% Perform EKI
my_eki500 = EKI(my_inverse,my_darcy,500);
my_eki500 = my_eki500.run_t(5);

my_eki1000 = EKI(my_inverse,my_darcy,1000);
my_eki1000 = my_eki1000.run_t(5);

my_eki5000 = EKI(my_inverse,my_darcy,5000);
my_eki5000 = my_eki5000.run_t(5);

plot_LMAP_vs_EKI(true_RTMflow,my_forward_mesh,my_inverse_mesh,my_darcy, ...
    lmap_mean,lmap_var, ...
    eki_mean500,eki_var500, ...
    eki_mean1000,eki_var1000, ...
    eki_mean5000,eki_var5000)

%% Perform MCMC
pool = gcp();
numWorkers = pool.NumWorkers;
my_mcmc = MCMC(my_inverse,my_darcy,100,numWorkers);
my_mcmc = my_mcmc.run_t(5);
