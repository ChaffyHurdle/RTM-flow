function obj = plotter(obj,u,C,h)

figure(1)

% Truth
subplot(2,2,1)
c_min = min(min(obj.inverse_class.u_true),min(u));
c_max = max(max(obj.inverse_class.u_true),max(u));
pdeplot(obj.inverse_class.fwd_mesh.nodes',...
        obj.inverse_class.fwd_mesh.elements', ...
        XYData=obj.inverse_class.u_true, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([c_min,c_max])
title("$u^\dagger$",'Interpreter','latex')
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off

% u_k
subplot(2,2,2)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',XYData=u, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([c_min,c_max])
title("$u_{map}$",'interpreter','latex')
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off
subplot(2,2,3)
pdeplot(obj.mesh_class.nodes',...
        obj.mesh_class.elements', ...
        XYData=h, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off

% h_k
title("$h$",'Interpreter','latex')
subplot(2,2,4)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',XYData=diag(C), ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([0,0.25])
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off
title("$C_{map}$",'interpreter','latex')

% figure(2)
% for i = 1:length(obj.physics_class.sensor_locs)
%     subplot(sqrt(length(obj.physics_class.sensor_locs)),sqrt(length(obj.physics_class.sensor_locs)),i)
%     plot(obj.inverse_class.data(i,:),'ko')
%     hold on
%     plot(candidate_RTM.pressure_data(i,:),'r*')
%     plot(G_u0(i,:),'b*')
%     hold off
%     ylim([obj.physics_class.p_0-0.1,obj.physics_class.p_I+0.1])
%     title("$D$ (black), G($u_0$) (blue), G($u_{map}$) (red)",'interpreter','latex')
% end

% posterior_samples = mvnrnd(u,C_post,25);
% figure(3)
% for i = 1:25
%     subplot(5,5,i)
%     pdeplot(obj.mesh_class.nodes',...
%             obj.mesh_class.elements', ...
%             XYData=posterior_samples(i,:), ...
%             XYStyle='interp',ColorMap="jet",Mesh="off")
%     clim([c_min,c_max])
% end

% Fog of war plot
figure(2)
centroid_x = obj.inverse_class.inv_mesh.centroids(:,1);
centroid_y = obj.inverse_class.inv_mesh.centroids(:,2);

% Compute distance-based transparency
% Initialize transparency as ones (fully opaque)
alpha_map = diag(C)/obj.inverse_class.matern_var;

% Create an alpha overlay
F = scatteredInterpolant(centroid_x, centroid_y, u', 'linear', 'none');
[xq, yq] = meshgrid(linspace(0,1,1000), ...
                    linspace(0,1,1000));
zq = F(xq, yq);

% Create an alpha overlay
F = scatteredInterpolant(centroid_x, centroid_y, alpha_map, 'linear', 'none');
[xq, yq] = meshgrid(linspace(0,1,1000), ...
                    linspace(0,1,1000));
alpha_overlay = F(xq, yq);

% Plot the alpha overlay as an image
im = imagesc(linspace(0,1,1000), linspace(0,1,1000), zq);
colormap jet
set(im, 'AlphaData', 1 - alpha_overlay); % Transparency decreases with alpha_map
set(im, 'AlphaDataMapping', 'none'); % Prevent scaling of alpha
clim([c_min,c_max])
colorbar
set(gca,'YDir','normal')
xlim([-1,1])
ylim([-1,1])
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off
drawnow

end
