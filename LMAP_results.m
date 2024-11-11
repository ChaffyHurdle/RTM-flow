addpath('Meshes')
meshes = {'p_ref.mat', 'e_ref.mat', 't_ref.mat', ...
    'p_new.mat', 'e_new.mat', 't_new.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end
%% Load

load("Results/u_iterations.mat","u_iterations");
load("Results/J_iterations.mat","J_iterations");
load("Results/scaled_data_misfit.mat","scaled_data_misfit");
load("Results/execution_times.mat","execution_times");
load("Results/C_map.mat","C_map");
load("Results/u_true.mat","u_true");

%% Mesh set up (avoiding inverse crimes)
my_forward_mesh = DelaunayMesh(p_ref,e_ref,t_ref);
my_inverse_mesh = DelaunayMesh(p_ref,e_ref,t_ref);
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~86% filled.
observation_times = linspace(0.1,0.9,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.92^2*mu*phi/(2*(p_I-p_0));

% Set N sensor locs (equally space)
sqrtN = 6;
sensor_locs_x = 1/(2*sqrtN) + linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = sensor_locs_x;
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

% Define true permeability and place within physics class
K_true = exp(u_true);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true, sensor_locs, observation_times,T);
my_pressure = Pressure(my_forward_mesh,my_darcy);

%% RTM set up (fine forward mesh)
true_RTMflow = RTMFlow(my_forward_mesh,my_darcy,my_pressure);
true_RTMflow = true_RTMflow.run();

t_index = find(true_RTMflow.times > my_darcy.observation_times(end),1)-1;

%% Plots
max_u = max(max(u_true,[],'all'),max(u_iterations(:,end),[],'all'));
min_u = min(min(u_true,[],'all'),min(u_iterations(:,end),[],'all'));

figure(1)
subplot(4,6,[1,2,7,8])
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements',XYData = u_true, ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
for j = 1:length(true_RTMflow.edge_data{t_index})
    plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
         [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
end
hold off
clim([min_u,max_u])
title('$u^{\dagger}$','interpreter','latex')

subplot(4,6,[3,4,9,10])
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements',XYData = u_iterations(:,end), ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
for j = 1:length(true_RTMflow.edge_data{t_index})
    plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
         [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
end
hold off
clim([min_u,max_u])
title("$u_{MAP}$, " + num2str(round(sum(execution_times),1)) + "s",'interpreter','latex')

subplot(4,6,[5,6,11,12])
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements',XYData = diag(C_map), ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
for j = 1:length(true_RTMflow.edge_data{t_index})
    plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
         [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
end
hold off
clim([0,0.25])
title("$C_{MAP}$",'interpreter','latex')

posterior_samples = mvnrnd(u_iterations(:,end),C_map,12);
for i = 1:12
    subplot(4,6,12+i)
    pdeplot(my_inverse_mesh.nodes',...
            my_inverse_mesh.elements', ...
            XYData=posterior_samples(i,:), ...
            XYStyle='interp',ColorMap="jet",Mesh="off")
    hold on
    for j = 1:length(true_RTMflow.edge_data{t_index})
        plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
             [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
    end
    hold off
    clim([min_u,max_u])
    title("Posterior sample")
end


figure(2)
subplot(1,2,1)
x = linspace(0,length(J_iterations)-1,length(J_iterations));
Y = [0.5*scaled_data_misfit; J_iterations-0.5*scaled_data_misfit]';
area(x,Y)
xlabel('Iteration')
ylabel('Cost function')
legend({'Data misfit contribution','Prior contribution'})
title("Cost function over iterations")

subplot(1,2,2)
plot(linspace(0,2*sqrtN^2*length(observation_times),1000),...
    chi2pdf(linspace(0,2*sqrtN^2*length(observation_times),1000),sqrtN^2*length(observation_times)))
hold on
xline(scaled_data_misfit(end))
hold off
title("Data misfit statistic vs. chi-squared distribution")


posterior_samples = mvnrnd(u_iterations(:,end),C_map,1000);
figure(3)
subplot(1,2,1)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements',XYData = u_true(:,end)>0.5 | u_true(:,end)<-0.5, ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
for j = 1:length(true_RTMflow.edge_data{t_index})
    plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
         [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
end
hold off

subplot(1,2,2)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements',XYData = mean(posterior_samples>0.5 | posterior_samples<-0.5), ...
    XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
for j = 1:length(true_RTMflow.edge_data{t_index})
    plot([my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),1), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),1)],...
         [my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,2),2), my_forward_mesh.nodes(true_RTMflow.edge_data{t_index}(j,3),2)],'w-','LineWidth',2)
end
hold off
