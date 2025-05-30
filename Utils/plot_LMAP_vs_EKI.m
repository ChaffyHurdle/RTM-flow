function p = plot_LMAP_vs_EKI(true_RTMflow,my_forward_mesh,my_inverse_mesh,my_darcy,my_lmap_seq,my_eki)

t_index = 5;

figure(5)
subplot(2,3,1)
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability), XYStyle='interp', ...
    ColorMap="turbo",Mesh="off")
clim([-1.5,1.5])
axis square;
title('$u^{\dagger}$','interpreter','latex')

subplot(2,3,2)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_lmap_seq.umap_seq(:,end),XYStyle='interp', ...
            ColorMap="turbo",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
hold off
clim([-1.5,1.5])
axis square;
title('LMAP mean')

subplot(2,3,3)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_lmap_seq.Cmap_seq(:,:,end)),XYStyle='interp', ...
            ColorMap="turbo",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
hold off
clim([0,max(my_inverse.C0_inv,[],'all')])
axis square;
title('LMAP variance')

subplot(2,3,4)
pdeplot(my_forward_mesh.nodes',my_forward_mesh.elements', ...
    XYData = log(my_darcy.permeability)*0, XYStyle='interp', ...
    ColorMap="turbo",Mesh="on")
hold on
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'wo','filled')
scatter(my_darcy.sensor_locs(:,1),my_darcy.sensor_locs(:,2),'ko')
hold off
clim([-1.5,1.5])
axis square;
axis on;
title('$Sensor locations$','interpreter','latex')


subplot(2,3,5)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = my_eki.ueki,XYStyle='interp', ...
            ColorMap="turbo",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
hold off
clim([-1.5,1.5])
axis square;
title('EKI mean')

subplot(2,3,6)
pdeplot(my_inverse_mesh.nodes',my_inverse_mesh.elements', ...
            XYData = diag(my_eki.Ceki),XYStyle='interp', ...
            ColorMap="turbo",Mesh="off")
hold on
plot_front(true_RTMflow,t_index)
hold off
clim([0,max(my_inverse.C0_inv,[],'all')])
axis square;
title('EKI variance')

end