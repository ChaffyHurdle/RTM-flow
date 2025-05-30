function im = plot_fog(my_lmap_seq,RTMflow_class,t)

x = linspace(0,1,1000);
y = x;
centroid_x = my_lmap_seq.inverse_class.inv_mesh.centroids(:,1);
centroid_y = my_lmap_seq.inverse_class.inv_mesh.centroids(:,2);
c_min = min(min(my_lmap_seq.inverse_class.u_true));
c_max = max(max(my_lmap_seq.inverse_class.u_true));
t_ind = find(RTMflow_class.times > my_lmap_seq.physics_class.observation_times(t),1,'first');

% Compute distance-based transparency
alpha_map = diag(my_lmap_seq.Cmap_seq(:,:,t))/my_lmap_seq.inverse_class.matern_var;

% Create an alpha overlay
F = scatteredInterpolant(centroid_x, centroid_y, my_lmap_seq.umap_seq(:,t), 'linear', 'none');
[xq, yq] = meshgrid(x,y);
zq = F(xq, yq);

% Create an alpha overlay
F = scatteredInterpolant(centroid_x, centroid_y, alpha_map, 'linear', 'none');
alpha_overlay = F(xq, yq);

% Plot the alpha overlay as an image
im = imagesc(x,y,zq);
colormap turbo
set(im, 'AlphaData', 1 - alpha_overlay); % Transparency decreases with alpha_map
set(im, 'AlphaDataMapping', 'none'); % Prevent scaling of alpha
%clim([c_min,c_max])
%colorbar
set(gca,'YDir','normal')
xlim([0,1])
ylim([0,1])
hold on
plot_front(RTMflow_class,t_ind)
hold off

end