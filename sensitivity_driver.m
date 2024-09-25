%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
addpath('Meshes')
meshes = {'p_new.mat', 'e_new.mat', 't_new.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end

%% Mesh set up
my_mesh = DelaunayMesh(p_new,e_new,t_new);

var_matern = 0.25; length_scale = 0.1; nu_matern = 1.5;
matern_args = [var_matern,length_scale,nu_matern];
my_inverse = Inversion(my_mesh,my_mesh,matern_args);

%% Generate u and h
my_inverse = my_inverse.generate_u();
h = my_inverse.generate_h();

%% Darcy setup
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~95% filled.
observation_times = linspace(0.2,0.8,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.95^2*mu*phi/(2*(p_I-p_0));

% Set sensor locs
sensor_locs_x = [0.25,0.5,0.75];
sensor_locs_y = [0.25,0.5,0.75];
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

% Define true permeability and place within physics class
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true, sensor_locs, observation_times,T);

%% Sensitivity setup
my_sensitivity = Sensitivity(my_mesh,my_darcy,my_inverse.u_true,h/4);
my_sensitivity = my_sensitivity.fwd_solves();
my_sensitivity = my_sensitivity.run();

%%
% i = 5;
% figure(1)
% pdeplot(my_mesh.nodes',my_mesh.elements',XYData= (my_sensitivity.pressures_u_plus_h(:,i)-my_sensitivity.pressures_u(:,i)) )
% hold on
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),2),'w.')
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),2),'wo')
% hold off
% figure(2)
% pdeplot(my_mesh.nodes',my_mesh.elements',XYData= my_sensitivity.is_moving_boundary_ob_times(:,i) )
% figure(3)
% pdeplot(my_mesh.nodes',my_mesh.elements',XYData= (my_sensitivity.pressures_u_plus_h(:,i)-my_sensitivity.pressures_u(:,i)).*my_sensitivity.is_moving_boundary_ob_times(:,i) )
% hold on
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),2),'w.')
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),2),'wo')
% hold off
% figure(4)
% pdeplot(my_mesh.nodes',my_mesh.elements',XYData= my_sensitivity.p_tildes(:,i) )
% hold on
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times(:,i)),2),'w.')
% plot(my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),1),my_mesh.nodes(boolean(my_sensitivity.is_moving_boundary_ob_times_u_plus_h(:,i)),2),'wo')
% hold off
% figure(5)
% pdeplot(my_mesh.nodes',my_mesh.elements',XYData= my_sensitivity.p_tildes(:,i).*my_sensitivity.is_moving_boundary_ob_times(:,i) )
% 
% 

% error_h = [];
% for i = 1:5
%     my_sensitivity = Sensitivity(my_mesh,my_darcy,my_inverse.u_true,h/(2^i));
%     my_sensitivity = my_sensitivity.fwd_solves();
%     my_sensitivity = my_sensitivity.run();
%     max_h = max(abs(my_sensitivity.h));
%     error_h = [error_h; [max_h, mean( (my_sensitivity.RTMflow_u_plus_h.pressure_data - my_sensitivity.RTMflow_u.pressure_data - my_sensitivity.p_tildes_at_sensors)/max_h,'all')]]
% end


