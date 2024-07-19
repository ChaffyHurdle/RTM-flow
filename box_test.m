%% Forward problem set up
my_mesh = DelaunayMesh(p_ref,e_ref,t_ref);

%% Inverse problem set up
matern_args = [0.25,0.1,1.5];
my_inverse = Inversion(my_mesh,matern_args);
my_inverse = my_inverse.generate_u();

%% Physics and pressure set up
mu = 0.1; phi = 1; thickness = 1; p_I = 6.0e5; p_0 = 1.0e5;
K_true = exp(my_inverse.u_meshcenters);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true);
my_pressure = Pressure(my_mesh,my_darcy);

%% Data extraction
% Times chosen so that in uniform permeability case, they are equally
% spaced between [0.1,0.9] (inclusive).
observation_times = linspace(0.1,0.9,7).^2*mu*phi/(2*(p_I-p_0));
sensor_locs_x = [0.2,0.4,0.6,0.8];
sensor_locs_y = [0.2,0.4,0.6,0.8];
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

%% RTM set up
my_RTMflow = RTMFlow(my_mesh,my_darcy,my_pressure,observation_times,sensor_locs);
my_RTMflow.visualise_class.is_plotting_volume = false;

figure(2)
pdeplot(my_inverse.DelaunayMesh.nodes',...
        my_inverse.DelaunayMesh.elements', ...
        XYData=my_inverse.u_meshcenters,XYStyle="flat",ColorMap="jet",Mesh="off")
colormap jet;
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
tic
my_RTMflow = my_RTMflow.run();
toc

%% Perform LMAP
my_inverse = my_inverse.generate_data(my_RTMflow.pressure_data,0.01);
my_lmap = LMAP(my_RTMflow, my_inverse);

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
% for i = 1:length(my_RTMflow.pressures)
%     figure(1)
%     pdeplot(my_RTMflow.Delaunay_mesh_class.nodes',...
%                 my_RTMflow.Delaunay_mesh_class.elements', ...
%                 XYData=my_RTMflow.pressures(:,i),ColorMap="jet",Mesh="on")
%     hold on
%     scatter(my_RTMflow.sensor_locs_on_mesh(:,1),my_RTMflow.sensor_locs_on_mesh(:,2),'wo','filled')
%     title("time elapsed: " + num2str(my_RTMflow.times(i)))
%     hold off
% end
% for i = 1:length(my_RTMflow.times)
%     figure(3)
%     pdeplot(my_RTMflow.Delaunay_mesh_class.nodes',...
%                 my_RTMflow.Delaunay_mesh_class.elements', ...
%                 XYData=my_RTMflow.filling_factors(:,i),ColorMap="jet",Mesh="on")
%     hold on
%     scatter(my_inverse.sensor_locs_on_mesh(:,1),my_inverse.sensor_locs_on_mesh(:,2),'wo','filled')
%     title("time elapsed: " + num2str(my_RTMflow.times(i)))
%     hold off
% end


