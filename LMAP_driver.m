%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
addpath('Meshes')
meshes = {'p_ref.mat', 'e_ref.mat', 't_ref.mat', ...
    'p_new.mat', 'e_new.mat', 't_new.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end
parpool('Threads', 10);

%% Mesh set up (avoiding inverse crimes)
my_forward_mesh = DelaunayMesh(p_ref,e_ref,t_ref);
my_inverse_mesh = DelaunayMesh(p_new,e_new,t_new);

%% Inverse problem set up
var_matern = 0.25; length_scale = 0.1; nu_matern = 1.5;
matern_args = [var_matern,length_scale,nu_matern];
my_inverse = Inversion(my_forward_mesh,my_inverse_mesh,matern_args);

% Generate true permeability to be recovered
my_inverse = my_inverse.generate_u();
my_inverse.u_true = my_inverse.u_true*0;
my_inverse.u_true(vecnorm( (my_forward_mesh.centroids-[0.25,0.3])' )<0.36^2) = 0.5;
my_inverse.u_true(vecnorm( (my_forward_mesh.centroids-[0.65,0.7])' )<0.36^2) = -0.6;
my_inverse.u_true(my_forward_mesh.centroids(:,2)>0.975) = 1;
my_inverse.u_true(my_forward_mesh.centroids(:,2)<0.025) = 1;
my_inverse.plot_u_true();

%% Physics and pressure set up
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~86% filled.
observation_times = linspace(0.15,0.85,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.86^2*mu*phi/(2*(p_I-p_0));

% Set N sensor locs (equally space)
sqrtN = 6;
sensor_locs_x = 1/(2*sqrtN) + linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = sensor_locs_x;
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
true_RTMflow = true_RTMflow.run();

%% Generate random data for inverse problem
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,0.005);

%% Perform LMAP
my_lmap = LMAP(my_inverse,my_darcy);
my_lmap = my_lmap.run();

u_iterations = my_lmap.u_iterations;
J_iterations = my_lmap.J_iterations;
execution_times = my_lmap.execution_times;
C_map = my_lmap.C_map;
u_true = my_inverse.u_true;
save("Results/u_iterations.mat","u_iterations");
save("Results/J_iterations.mat","J_iterations");
save("Results/execution_times.mat","execution_times");
save("Results/C_map.mat","C_map");
save("Results/u_true.mat","u_true");


% %% To be tidied
% subplot(1,3,1)
% max_c = max(max((1+my_lmap.alpha)*my_lmap.tildePmat,[],'all'),max(my_lmap.tildePmat2,[],'all'));
% imagesc((1+my_lmap.alpha)*my_lmap.tildePmat)
% colorbar
% clim([0,max_c])
% subplot(1,3,2)
% imagesc(my_lmap.tildePmat2)
% colorbar
% clim([0,max_c])
% subplot(1,3,3)
% imagesc((1+my_lmap.alpha)*my_lmap.tildePmat-my_lmap.tildePmat2)
% colorbar
% clim([0,max_c])
% 
% figure(1)
% for i = 0:4
%     for j = 1:9
%         subplot(9,5,9*i+j)
%         pdeplot(my_inverse_mesh.nodes',...
%         my_inverse_mesh.elements', ...
%         XYData=my_lmap.R(:,9*i+j),XYStyle='interp',ColorMap="jet",Mesh="off")
%         hold on
%         scatter(my_lmap.physics_class.sensor_locs(:,1),my_lmap.physics_class.sensor_locs(:,2),'wo','filled')
%         scatter(my_lmap.physics_class.sensor_locs(my_lmap.i_vec(9*i+j),1),my_lmap.physics_class.sensor_locs(my_lmap.i_vec(9*i+j),2),'ko')
%         hold off
%         %caxis([min(min(my_lmap.representers)),max(max(my_lmap.representers))])
%     end
%     subplot(9,5,i+1)
%     title("Time: " + num2str(i+1))
% end
% figure(2)
% for i = 0:4
%     for j = 1:9
%         subplot(9,5,9*i+j)
%         pdeplot(my_inverse_mesh.nodes',...
%         my_inverse_mesh.elements', ...
%         XYData=my_lmap.Q(:,9*i+j),XYStyle='interp',ColorMap="jet",Mesh="off")
%         %caxis([min(min(my_lmap.representers)),max(max(my_lmap.representers))])
%         hold on
%         scatter(my_lmap.physics_class.sensor_locs(:,1),my_lmap.physics_class.sensor_locs(:,2),'wo','filled')
%         scatter(my_lmap.physics_class.sensor_locs(my_lmap.i_vec(9*i+j),1),my_lmap.physics_class.sensor_locs(my_lmap.i_vec(9*i+j),2),'ko')
%         hold off
%     end
%     subplot(9,5,i+1)
%     title("Time: " + num2str(i+1))
% end
% 
% k = 3;
% disp(num2str(my_lmap.physics_class.observation_times(my_lmap.j_vec(k))))
% for i = 1:length(my_lmap.RTMflow_class.times)
%     figure(1)
%     pdeplot(my_inverse_mesh.nodes',...
%                my_inverse_mesh.elements', ...
%                XYData=my_lmap.p_tildes(:,i,k),XYStyle='interp',ColorMap="jet",Mesh="off")
%     clim([min(my_lmap.p_tildes(:,:,k),[],'all'),max(my_lmap.p_tildes(:,:,k),[],'all')])
%     hold on
%     scatter(my_lmap.physics_class.sensor_locs(:,1),my_lmap.physics_class.sensor_locs(:,2),'wo','filled')
%     scatter(my_lmap.physics_class.sensor_locs(my_lmap.i_vec(k),1),my_lmap.physics_class.sensor_locs(my_lmap.i_vec(k),2),'ko')
%     hold off
%     title(num2str(my_lmap.RTMflow_class.times(i)))
% end
% 
% for i = 1:length(my_lmap.RTMflow_class.times)
%     figure(1)
%     pdeplot(my_inverse_mesh.nodes',...
%                my_inverse_mesh.elements', ...
%                XYData=my_lmap.p_tilde0(:,i),XYStyle='interp',ColorMap="jet",Mesh="off")
%     clim([0,1])
%     title(num2str(my_lmap.RTMflow_class.times(i)))
% end