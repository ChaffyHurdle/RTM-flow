%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
parpool('Threads');
addpath('Utils')

%% Set case study number and load meshes and sensor locs
case_study = 1;

% Add all relevant meshes and sensor locations
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

% Load inlet/outlet function.
% polar = true for radial flow patterns (curvature computed differently)
[inlet_func,vent_func,polar] = load_inlets_vents(case_study);

%% Mesh set up (avoiding inverse crimes)
my_forward_mesh = DelaunayMesh(p_fwd,e_fwd,t_fwd);
my_inverse_mesh = DelaunayMesh(p_inv,e_inv,t_inv);

%% Load parameters (prior and Darcy flow). These can also be set manually.
[matern_args,mu,phi,thickness,p_I,p_0,observation_times,T] = load_parameters(case_study);

%% Inverse problem set up (~20 secs)
my_inverse = Inversion(my_forward_mesh,my_inverse_mesh,matern_args);

%% Generate new true permeability, if wanted

% Pull previous example:
% u_mat = load("path_to_your_log_permeability_file.mat");
% my_inverse.u_true = u_mat.u;

% Generate new permeability
while true
    my_inverse = my_inverse.generate_u();
    my_inverse.plot_u_true();
    key_input = input("Confirm ('y') or sample again (any other key)");
    if key_input == 'y'
        break
    end
end

%% Physics and pressure set up
sensor_locs = sensor_locs_25; % set to personal preference
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, ...
    p_I, p_0, inlet_func, vent_func, ...
    K_true, sensor_locs, observation_times,T);
my_pressure = Pressure(my_forward_mesh,my_darcy);

% Plot experimental setup
figure(1)
my_inverse.plot_u_true()
hold on
scatter(sensor_locs(:,1),sensor_locs(:,2),'white','filled')
xlim([min(my_forward_mesh.nodes(:,1)),max(my_forward_mesh.nodes(:,1))])
ylim([min(my_forward_mesh.nodes(:,2)),max(my_forward_mesh.nodes(:,2))])
hold off

%% RTM true simulation (on fine forward mesh)
true_RTMflow = RTMFlow(my_forward_mesh,my_darcy,my_pressure,polar);
true_RTMflow = true_RTMflow.run(inf);

if T > true_RTMflow.time
    disp('T greater than tau')
end
%% Generate random data for inverse problem
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,0.005);

%% Perform LMAP (all times)
my_lmap = LMAP(my_inverse,my_darcy,1e3,2,0.03,0.03,polar);
my_lmap = my_lmap.run();

save(strcat(folder_path,'/lmap_seq.mat'), 'my_lmap')

plot_LMAP_seq(my_forward_mesh, my_inverse_mesh, my_darcy, ...
    true_RTMflow,my_lmap)

%% Plot push forward example

u_samples = mvnrnd(my_lmap.umap_seq(:,end),my_lmap.Cmap_seq(:,:,end),1000);
[pressures,flow_fronts] = push_forward(u_samples, my_darcy, my_inverse_mesh);
perturbed_pressures = pressures + normrnd(0,1,size(pressures)).*sqrt(my_inverse.Sigma(:)');
plot_push_forward(perturbed_pressures, flow_fronts, my_darcy, my_inverse_mesh, true_RTMflow)

%% Perform EKI
my_eki = EKI(my_inverse,my_darcy,1000);
my_eki = my_eki.run_t(5);

plot_LMAP_vs_EKI(true_RTMflow,my_forward_mesh,my_inverse_mesh,my_darcy,my_lmap,my_eki)

%% Perform MCMC
pool = gcp();
numWorkers = pool.NumWorkers;
my_mcmc = MCMC(my_inverse,my_darcy,100,numWorkers);
my_mcmc = my_mcmc.run_t(5);

%%
% figure(3)
% for j = round(14*length(true_RTMflow.edge_data)/15):length(true_RTMflow.edge_data)
%     plot(0,0)
%     hold on
%     for k = 1:length(true_RTMflow.edge_data{j})
%         plot([true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1), ...
%             true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1)],...
%             [true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2), ...
%             true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2)],'k-')
%         % plot([(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1))/2, ...
%         %       (true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1))/2 + true_RTMflow.edge_data{j}(k,4)/20], ...
%         %      [(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2))/2, ...
%         %       (true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2))/2 + true_RTMflow.edge_data{j}(k,5)/20],'k-')
%     end
%     scatter(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(:,2),1),...
%             true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(:,2),2),'b*')
%     % scatter(candidate_RTM.Delaunay_mesh_class.nodes(candidate_RTM.moving_boundary(:,j),1),...
%     %         candidate_RTM.Delaunay_mesh_class.nodes(candidate_RTM.moving_boundary(:,j),2),'bo')
%     tester = (boolean(true_RTMflow.active_nodes(:,j)) & boolean(true_RTMflow.Dirichlet_nodes(:,j))) | boolean(true_RTMflow.pressure_class.is_vent);
%     scatter(true_RTMflow.Delaunay_mesh_class.nodes(tester,1),...
%             true_RTMflow.Delaunay_mesh_class.nodes(tester,2),'bo')
%     plot(cos(0:pi/50:2*pi), sin(0:pi/50:2*pi),'r-');
%     xlim([-1.1,1.1])
%     ylim([-1.1,1.1])
%     hold off
%     drawnow
% end
