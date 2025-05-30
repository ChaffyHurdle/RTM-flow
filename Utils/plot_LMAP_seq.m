function c1 = plot_LMAP_seq(my_forward_mesh, my_inverse_mesh, my_darcy, ...
    true_RTMflow,my_lmap)

% Create tiled layout
t = tiledlayout(2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

% Preallocate axes handles
axTop = gobjects(1, 6);
axBottom = gobjects(1, 6);

% First plot
axTop(1) = nexttile(1);
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability), XYStyle='interp', ...
    Mesh="off",ColorBar='off')
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
axis square;
title('$u^{\dagger}$','interpreter','latex')

% Plot first row (2 to 6)
for i = 2:6
    axTop(i) = nexttile(i);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;

    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_lmap.umap_seq(:,i-1),XYStyle='interp',Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index)
    scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
    scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
    hold off
    axis square;
    title(sprintf('$\\mu_{%d}$', i), 'Interpreter', 'latex');
    axis off;
    colormap('turbo');
    clim([-1.5 1.5]);
end

% Plot second row (7 to 12)
for i = 2:6
    axBottom(i) = nexttile(i + 6);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;
    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_lmap.Cmap_seq(:,:,i-1)),XYStyle='interp', ...
            Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index)
    scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
    scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
    hold off
    axis square;
    axis off;
    colormap('turbo');
    clim([0 0.25]);
end

% Add shared colorbar for top row (using the last axis of the row)
cb1 = colorbar(axTop(6), 'Location', 'eastoutside');  % Colorbar for top row

% Add shared colorbar for bottom row (using the last axis of the row)
cb2 = colorbar(axBottom(6), 'Location', 'eastoutside');  % Colorbar for bottom row

% Now manually adjust the colorbar positions to match row heights
% Get the position of the last axes in each row (for height alignment)
topPosition = axTop(6).Position;  % Get position of the last axis in top row
bottomPosition = axBottom(6).Position;  % Get position of the last axis in bottom row

% Adjust colorbar 1 (for top row)
cb1.Position(4) = topPosition(4);  % Match height with the top row
cb1.Position(2) = topPosition(2);  % Align the vertical position

% Adjust colorbar 2 (for bottom row)
cb2.Position(4) = bottomPosition(4);  % Match height with the bottom row
cb2.Position(2) = bottomPosition(2);  % Align the vertical position

cb1.Limits = [-1.5 1.5];
cb2.Limits = [0 0.25];

end

