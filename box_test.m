%% clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
addpath('Meshes')
meshes = {'p_ref.mat', 'e_ref.mat', 't_ref.mat', ...
    'p.mat', 'e.mat', 't.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end
parpool('local', 10);

%% Mesh set up (avoid inverse crimes)
my_forward_mesh = DelaunayMesh(p_ref,e_ref,t_ref);
my_inverse_mesh = DelaunayMesh(p,e,t);

%% Inverse problem set up
var_matern = 0.25; length_scale = 0.1; nu_matern = 1.5;
matern_args = [var_matern,length_scale,nu_matern];
my_inverse = Inversion(my_forward_mesh,my_inverse_mesh,matern_args);
my_inverse = my_inverse.generate_u();
my_inverse.plot_u_true();

%% Physics and pressure set up
mu = 0.1; phi = 1; thickness = 1; p_I = 6.0e5; p_0 = 1.0e5;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~95% filled.
observation_times = linspace(0.1,0.9,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.95^2*mu*phi/(2*(p_I-p_0));

% Set sensor locs.
sensor_locs_x = [0.25,0.5,0.75];
sensor_locs_y = [0.25,0.5,0.75];
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

% Define true permeability and place within physics class
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true, sensor_locs, observation_times,T);
my_pressure = Pressure(my_forward_mesh,my_darcy);

%% RTM set up (fine forward mesh)
true_RTMflow = RTMFlow(my_forward_mesh,my_darcy,my_pressure);
% true_RTMflow.visualise_class.is_plotting_volume = true;
true_RTMflow = true_RTMflow.run();
% true_RTMflow.flow_plots();

%% Generate data for inverse problem
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,0.01);

%% Perform LMAP
my_lmap = LMAP(my_inverse,my_darcy);
my_lmap = my_lmap.run();

figure(1)
for i = 0:4
    for j = 1:9
        subplot(9,5,9*i+j)
        pdeplot(my_inverse_mesh.nodes',...
        my_inverse_mesh.elements', ...
        XYData=my_lmap.representers(:,9*i+j),XYStyle='interp',ColorMap="jet",Mesh="off")
        hold on
        scatter(my_lmap.RTMflow_class.sensor_locs_on_mesh(:,1),my_lmap.RTMflow_class.sensor_locs_on_mesh(:,2),'wo','filled')
        hold off
        %caxis([min(min(my_lmap.representers)),max(max(my_lmap.representers))])
    end
end
figure(2)
for i = 0:4
    for j = 1:9
        smoothed_reps = my_inverse.C0_inv * (my_lmap.representers .* my_inverse_mesh.element_areas);
        subplot(9,5,9*i+j)
        pdeplot(my_inverse_mesh.nodes',...
        my_inverse_mesh.elements', ...
        XYData=smoothed_reps(:,9*i+j),XYStyle='interp',ColorMap="jet",Mesh="off")
        %caxis([min(min(my_lmap.representers)),max(max(my_lmap.representers))])
        hold on
        scatter(my_lmap.RTMflow_class.sensor_locs_on_mesh(:,1),my_lmap.RTMflow_class.sensor_locs_on_mesh(:,2),'wo','filled')
        hold off
    end
end


%% Execute solver
% parallel_pressure = cell(1,3);
% tic
% parfor i=1:3
%     disp(i)
%     visc = [0.09,0.1,0.11]
%     my_darcy = Physics(visc(i), 0.67, 1, @permeability);
%     my_RTMflow = RTMFlow(my_mesh,my_darcy,my_pressure,observation_times,sensor_locs);
%     my_RTMflow = my_RTMflow.run();
%     parallel_pressure{i} = my_RTMflow.pressure_data;
% end
% toc
k = 4;
for i = 1:width(my_lmap.lambdas(:,:,k))
    figure(1)
    pdeplot(my_inverse_mesh.nodes',...
                my_inverse_mesh.elements', ...
                XYData=my_lmap.lambdas(:,i,k),XYStyle='interp',ColorMap="jet",Mesh="off")
    hold on
    scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
    title("time elapsed: " + num2str(my_lmap.RTMflow_class.times(i)))
    hold off

%     figure(2)
%     pdeplot(my_inverse_mesh.nodes',...
%                 my_inverse_mesh.elements', 'flowdata',...
%                 my_lmap.grad_lambdas(:,:,i,k),XYStyle='interp',ColorMap="jet",Mesh="off")
%     hold on
%     scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
%     title("time elapsed: " + num2str(true_RTMflow.times(i)))
%     hold off
%     xlim([0,1]); ylim([0,1]);
end
