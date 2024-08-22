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
my_inverse.plot_u_true();
h = my_inverse.generate_h();

%% Darcy setup
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~95% filled.
observation_times = linspace(0.1,0.9,5).^2*mu*phi/(2*(p_I-p_0));
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
my_sensitivity = Sensitivity(my_mesh,my_darcy,my_inverse.u_true,h);
my_sensitivity = my_sensitivity.fwd_solves();

%% Plot differences in moving front and pressure

figure(1)
k = 1;
for i = my_darcy.observation_times
    subplot(1,5,k)
    u_RTMflow = my_sensitivity.RTMflow_u;
    [ d, ix ] = min( abs( u_RTMflow.times-i ) );
    front = boolean(u_RTMflow.active_nodes(:,ix) & u_RTMflow.Dirichlet_nodes(:,ix) & ~u_RTMflow.pressure_class.is_inlet);
    front_nodes = sortrows(my_mesh.nodes(front,:),2);
    plot(front_nodes(:,1),front_nodes(:,2))
    
    hold on
    u_plus_h_RTMflow = my_sensitivity.RTMflow_u_plus_h{1};
    [ d, ix ] = min( abs( u_plus_h_RTMflow.times-i ) );
    front = boolean(u_plus_h_RTMflow.active_nodes(:,ix) & u_plus_h_RTMflow.Dirichlet_nodes(:,ix) & ~u_plus_h_RTMflow.pressure_class.is_inlet);
    front_nodes = sortrows(my_mesh.nodes(front,:),2);
    plot(front_nodes(:,1),front_nodes(:,2))
    hold off

    xlim([0,1]);ylim([0,1])
    k = k+1;
end

figure(2)
u_RTMflow = my_sensitivity.RTMflow_u;
max_error = zeros(1,4);
k = 1;
for i=1:4
    u_plus_h_RTMflow = my_sensitivity.RTMflow_u_plus_h{i};
    for j=my_darcy.observation_times
        subplot(4,5,k)
        [ d, ix_u_plus_h ] = min( abs( u_plus_h_RTMflow.times-j ) );
        [ d, ix_u ] = min( abs( u_RTMflow.times-j ) );
        h_norm = max(2^(-i)*my_sensitivity.h);
        pdeplot(my_mesh.nodes',...
               my_mesh.elements', ...
               XYData=(u_plus_h_RTMflow.pressures(:,ix_u_plus_h) - u_RTMflow.pressures(:,ix_u))/h_norm,XYStyle='interp',ColorMap="jet",Mesh="off")
        colorbar; caxis([-1/2 1/2]);
        k = k + 1;
    end
    max_error(i) = (u_plus_h_RTMflow.pressures(100,ix_u_plus_h) - u_RTMflow.pressures(100,ix_u))/h_norm;
end


