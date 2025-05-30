%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
parpool('Threads');

%% Set case study number and load data
case_study = 3;
folder_path = strcat('Case',num2str(case_study));  % Change this
files = dir(fullfile(folder_path, '*.mat'));

all_data = struct();  % to hold everything
for k = 1:length(files)
    file_path = fullfile(folder_path, files(k).name);
    [~, name, ~] = fileparts(files(k).name);
    tmp = load(file_path);
    var_names = fieldnames(tmp);
    assignin('base', name, tmp.(var_names{1}));
end
[inlet_func,vent_func,polar] = load_inlets_vents(case_study);

%% Mesh set up (avoiding inverse crimes)
my_forward_mesh = DelaunayMesh(p_fwd,e_fwd,t_fwd);
my_inverse_mesh = DelaunayMesh(p_inv,e_inv,t_inv);

%% Inverse problem set up
var_matern = 0.25; length_scale = 0.1; nu_matern = 1.5;
matern_args = [var_matern,length_scale,nu_matern];
my_inverse = Inversion(my_forward_mesh,my_inverse_mesh,matern_args);

%% Generate new true permeability, if wanted
while true
    my_inverse = my_inverse.generate_u();
    my_inverse.plot_u_true();
    key_input = input("Confirm ('y') or sample again (any other key)");
    if key_input == 'y'
        break
    end
end

%% Physics and pressure set up
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~86% filled.
observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.92^2*mu*phi/(2*(p_I-p_0));

% Set N sensor locs (equally spaced)
sqrtN = 5;
sensor_locs_x = 1/(2*sqrtN) + linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = sensor_locs_x;
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];
save(strcat("Case1/sensor_locs_",num2str(sqrtN^2),".mat"),"sensor_locs")
disp([sqrtN^2,length(observation_times)])

n_per_turn = 16;
n_radii = 6;
turns = 2*pi * (0:(n_per_turn-1)) / n_per_turn;
radii = linspace(0.35,0.9,n_radii);
sensor_locs_radii = reshape(repmat(radii,n_per_turn,1),1,[]);
sensor_locs_theta = repmat(turns,1,n_radii) + reshape(repmat(double(~mod(1:n_radii,2))*(turns(2)-turns(1))/2,n_per_turn,1),1,[]);
sensor_locs_x = sensor_locs_radii .* cos(sensor_locs_theta);
sensor_locs_y = sensor_locs_radii .* sin(sensor_locs_theta);
sensor_locs = [sensor_locs_x' sensor_locs_y'];

sqrtN = 6;
sensor_locs_x = 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = 0.25 + 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

sensor_locs_x_prong = 0.5 + 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y_prong = sort(unique(sensor_locs_y));
for i = 1:sqrtN
    for j = 1:sqrtN/2
        %prong_upper_y = 0.5*sensor_locs_x_prong(i) + (sensor_locs_y_prong(j)-0.25);
        prong_lower_y = -0.5*sensor_locs_x_prong(i) + (sensor_locs_y_prong(j)+0.25);
        sensor_locs = [sensor_locs; sensor_locs_x_prong(i) 1-prong_lower_y];
        sensor_locs = [sensor_locs; sensor_locs_x_prong(i) prong_lower_y];
    end
end

% sensor_locs_x = 0.5 + 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
% sensor_locs_y = 0.5*sensor_locs_y;
% %save(strcat("Case1/sensor_locs_",num2str(sqrtN^2),".mat"),"sensor_locs")
% disp([sqrtN^2,length(observation_times)])

figure(1)
my_inverse.plot_u_true()
hold on
scatter(sensor_locs(:,1),sensor_locs(:,2),'red','filled')
hold off


% Define true permeability and place within physics class
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, ...
    p_I, p_0, inlet_func, vent_func, ...
    K_true, sensor_locs, observation_times,T);
my_pressure = Pressure(my_forward_mesh,my_darcy);

%% RTM true simulation (on fine forward mesh)
true_RTMflow = RTMFlow(my_forward_mesh,my_darcy,my_pressure,polar);
true_RTMflow = true_RTMflow.run(inf);

%% Generate random data for inverse problem
my_inverse = my_inverse.generate_data(true_RTMflow.pressure_data,0.005);

%% Perform LMAP (all times)
my_lmap = LMAP(my_inverse,my_darcy,1e3,2,0.03,0.03,polar);
my_lmap = my_lmap.run();

figure(4)
ax(1,1) = subplot(2,6,1);
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability), XYStyle='interp', ...
    Mesh="off",ColorBar='off')
%axis square;
axis off;
title('$u^{\dagger}$','interpreter','latex')

% Plot first row (2 to 6)
for i = 2:6
    figure(4)
    ax(1,i) = subplot(2,6,i);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;

    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_lmap.umap_seq(:,i-1),XYStyle='interp',Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index);
    hold off
    %axis square;
    axis off;
    title(sprintf('$\\bar{u}_{%d}$', i-1), 'Interpreter', 'latex');
    colormap('turbo');
    clim([-1.5 1.5]);
end
colormap(ax(1,1),'turbo');
pos = get(subplot(2,6,6),'Position');
h = colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  pos(3)/10  pos(4)]);


subplot(2,6,7)
colormap('gray');
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability)*0, XYStyle='interp', ...
    Mesh="on",ColorBar='off')
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
%axis square;
axis on;

% Plot second row (7 to 12)
for i = 2:6
    ax(2,i) = subplot(2,6,i+6);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;
    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_lmap.Cmap_seq(:,:,i-1)),XYStyle='interp', ...
            Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index);
    hold off
    %axis square;
    axis off;
    title(sprintf('$\\mathcal{C}_{%d}$', i-1), 'Interpreter', 'latex');
    colormap('turbo');
    clim([0 0.25]);
end
colormap(ax(2,6),'turbo');
pos = get(subplot(2,6,12),'Position');
h = colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  pos(3)/10  pos(4)]);




%% Plot simple example
plot_LMAP_seq(my_forward_mesh, my_inverse_mesh, my_darcy, ...
    true_RTMflow,my_lmap);

u_samples = mvnrnd(my_lmap.u_map,my_lmap.C_map,100);
[pressures,flow_fronts] = push_forward(u_samples, my_darcy, my_inverse_mesh);
perturbed_pressures = pressures + normrnd(0,1,size(pressures)).*sqrt(my_inverse.Sigma(:)');
plot_push_forward(perturbed_pressures, flow_fronts, my_darcy, my_inverse_mesh, true_RTMflow)


figure(3)
for j = round(14*length(true_RTMflow.edge_data)/15):length(true_RTMflow.edge_data)
    plot(0,0)
    hold on
    for k = 1:length(true_RTMflow.edge_data{j})
        plot([true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1), ...
            true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1)],...
            [true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2), ...
            true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2)],'k-')
        % plot([(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1))/2, ...
        %       (true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),1) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),1))/2 + true_RTMflow.edge_data{j}(k,4)/20], ...
        %      [(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2))/2, ...
        %       (true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,2),2) + true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(k,3),2))/2 + true_RTMflow.edge_data{j}(k,5)/20],'k-')
    end
    scatter(true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(:,2),1),...
            true_RTMflow.Delaunay_mesh_class.nodes(true_RTMflow.edge_data{j}(:,2),2),'b*')
    % scatter(candidate_RTM.Delaunay_mesh_class.nodes(candidate_RTM.moving_boundary(:,j),1),...
    %         candidate_RTM.Delaunay_mesh_class.nodes(candidate_RTM.moving_boundary(:,j),2),'bo')
    tester = (boolean(true_RTMflow.active_nodes(:,j)) & boolean(true_RTMflow.Dirichlet_nodes(:,j))) | boolean(true_RTMflow.pressure_class.is_vent);
    scatter(true_RTMflow.Delaunay_mesh_class.nodes(tester,1),...
            true_RTMflow.Delaunay_mesh_class.nodes(tester,2),'bo')
    plot(cos(0:pi/50:2*pi), sin(0:pi/50:2*pi),'r-');
    xlim([-1.1,1.1])
    ylim([-1.1,1.1])
    hold off
    drawnow
end

%% Perform EKI
my_eki = EKI(my_inverse,my_darcy);
my_eki = my_eki.run_t(5);

my_eki = EKI(my_inverse,my_darcy,100);
my_eki = my_eki.run_t(5);

%% Perform MCMC
pool = gcp();
numWorkers = pool.NumWorkers;
my_mcmc = MCMC(my_inverse,my_darcy,100,numWorkers);
my_mcmc = my_mcmc.run_t(5);

%% Save data
% u_iterations = my_lmap.u_iterations;
% J_iterations = my_lmap.J_iterations;
% scaled_data_misfit = my_lmap.scaled_data_misfit;
% execution_times = my_lmap.execution_times;
% C_map = my_lmap.C_map;
% u_true = my_inverse.u_true;
% save("Results/u_iterations.mat","u_iterations");
% save("Results/J_iterations.mat","J_iterations");
% save("Results/scaled_data_misfit.mat","scaled_data_misfit");
% save("Results/execution_times.mat","execution_times");
% save("Results/C_map.mat","C_map");
% save("Results/u_true.mat","u_true");

figure(5)
subplot(2,3,1)
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability), XYStyle='interp', ...
    ColorMap="jet",Mesh="off")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([-1.5,1.5])
axis square;
title('$u^{\dagger}$','interpreter','latex')

subplot(2,3,2)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_lmap.umap_seq(:,i),XYStyle='interp', ...
            ColorMap="jet",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([-1.5,1.5])
axis square;
title('LMAP mean')

subplot(2,3,3)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_lmap.Cmap_seq(:,:,i)),XYStyle='interp', ...
            ColorMap="jet",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([0,0.25])
axis square;
title('LMAP variance')

subplot(2,3,5)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_eki.ueki,XYStyle='interp', ...
            ColorMap="jet",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([-1.5,1.5])
axis square;
title('EKI mean')

subplot(2,3,6)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_eki.Ceki),XYStyle='interp', ...
            ColorMap="jet",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([0, 0.25])
axis square;
title('EKI variance')
