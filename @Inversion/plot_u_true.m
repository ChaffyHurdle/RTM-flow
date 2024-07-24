function obj = plot_u_true(obj)

figure(1)
clf;

pdeplot(obj.fwd_mesh.nodes',obj.fwd_mesh.elements', ...
    XYData=obj.u_true,XYStyle='flat',ColorMap='jet',Mesh='on')
axis equal
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
title("True log-permeability")

end