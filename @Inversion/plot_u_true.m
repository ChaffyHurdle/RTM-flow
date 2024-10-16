function obj = plot_u_true(obj)

figure(1)
clf;

pdeplot(obj.fwd_mesh.nodes',obj.fwd_mesh.elements', ...
    XYData=obj.u_true,XYStyle='interp',ColorMap='jet',Mesh='off')
axis equal
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
title("True log-permeability")
clim([-1.5,1.5])

end