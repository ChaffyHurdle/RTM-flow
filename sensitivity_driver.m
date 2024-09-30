%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
addpath('Meshes')
meshes = {'p_ref.mat', 'e_ref.mat', 't_ref.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end

%% Mesh set up
my_mesh = DelaunayMesh(p_ref,e_ref,t_ref);

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
my_sensitivity = Sensitivity(my_mesh,my_darcy,my_inverse.u_true,h/2);
my_sensitivity = my_sensitivity.fwd_solves();
my_sensitivity = my_sensitivity.run();
my_sensitivity = my_sensitivity.plot_sensitivity();

