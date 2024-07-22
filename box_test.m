%% clearing and loading
clear all
addpath('Meshes')
meshes = {'p_ref.mat', 'e_ref.mat', 't_ref.mat', ...
    'p.mat', 'e.mat', 't.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end
parpool('IdleTimeout', 120);

%% Forward problem set up
my_refined_mesh = DelaunayMesh(p,e,t);

%% Inverse problem set up
matern_args = [0.25,0.1,1.5];
my_inverse = Inversion(my_refined_mesh,matern_args);
my_inverse = my_inverse.generate_u();
my_inverse.plot_u_true();

%% Physics and pressure set up
mu = 0.1; phi = 1; thickness = 1; p_I = 6.0e5; p_0 = 1.0e5;
K_true = exp(my_inverse.u_meshcenters);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true);
my_pressure = Pressure(my_refined_mesh,my_darcy);

%% Data extraction
% Times chosen so that in uniform permeability case, they are equally
% spaced between [0.1,0.9] (inclusive).
observation_times = linspace(0.1,0.9,7).^2*mu*phi/(2*(p_I-p_0));
T = 0.95^2*mu*phi/(2*(p_I-p_0));

sensor_locs_x = [0.2,0.4,0.6,0.8];
sensor_locs_y = [0.2,0.4,0.6,0.8];
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

%% RTM set up
true_RTMflow = RTMFlow(my_refined_mesh,my_darcy, ...
    my_pressure,observation_times,sensor_locs,T);
true_RTMflow.visualise_class.is_plotting_volume = false;

true_RTMflow = true_RTMflow.run();
%true_RTMflow.flow_plots();

%% Perform LMAP
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,0.01);
my_lmap = LMAP(true_RTMflow, my_inverse);
my_lmap = my_lmap.compute_all_lambdas();

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
    pdeplot(true_RTMflow.Delaunay_mesh_class.nodes',...
                true_RTMflow.Delaunay_mesh_class.elements', ...
                XYData=my_lmap.lambdas(:,end-1,k),ColorMap="jet",Mesh="on")
    hold on
    scatter(true_RTMflow.sensor_locs_on_mesh(:,1),true_RTMflow.sensor_locs_on_mesh(:,2),'wo','filled')
    title("time elapsed: " + num2str(true_RTMflow.times(i)))
    hold off
end


