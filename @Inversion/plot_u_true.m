function obj = plot_u_true(obj)

figure(1)
clf;

pdeplot(obj.DelaunayMesh.nodes',obj.DelaunayMesh.elements', ...
    XYData=obj.u_meshcenters,XYStyle='flat',ColorMap='jet',Mesh='on')
axis equal
caxis([-1.75 1.75]);
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
title("True log-permeability")

end