function c1 = plot_LMAP_seq(my_forward_mesh, my_inverse_mesh, my_darcy, ...
    true_RTMflow,my_lmap)

c_min = min(log(my_darcy.permeability));
c_max = max(log(my_darcy.permeability));

%% True permeability
figure(4)
ax(1,1) = subplot(2,6,1);
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability), XYStyle='interp', ...
    Mesh="off",ColorBar='off')
clim([c_min c_max]);
axis off;
title('$u^{\dagger}$','interpreter','latex')

%% Mean estimates
for i = 2:6
    figure(4)
    ax(1,i) = subplot(2,6,i);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;
    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_lmap.umap_seq(:,i-1),XYStyle='interp',Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index);
    hold off
    axis off;
    title(sprintf('$\\bar{u}_{MAP}^{(%d)}$', i-1), 'Interpreter', 'latex');
    colormap('turbo');
    clim([c_min c_max]);
end
colormap(ax(1,1),'turbo');
pos = get(ax(1,6),'Position');
h = colorbar('Position', [pos(1)+(1+1/20)*pos(3)  pos(2)  pos(3)/10  pos(4)]);

%% Sensor configuration
ax(2,1) = subplot(2,6,7);
colormap('gray');
plot(1,1)
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
axis on;
xlim([0,1])
ylim([0,1])
title('Sensor locations','Interpreter','latex')

%% Variance estimates
for i = 2:6
    ax(2,i) = subplot(2,6,i+6);
    t_index = find(true_RTMflow.times > my_darcy.observation_times(i-1),1)-1;
    pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_lmap.Cmap_seq(:,:,i-1)),XYStyle='interp', ...
            Mesh="off",ColorBar='off')
    hold on
    plot_front(true_RTMflow,t_index);
    hold off
    axis off;
    title(sprintf('diag$(\\mathcal{C}_{MAP}^{(%d)})$', i-1), 'Interpreter', 'latex');
    colormap('turbo');
    clim([0 my_lmap.inverse_class.matern_var]);
end
colormap(ax(2,6),'turbo');
pos = get(ax(2,6),'Position');
h = colorbar('Position', [pos(1)+(1+1/20)*pos(3)  pos(2)  pos(3)/10  pos(4)]);

end

