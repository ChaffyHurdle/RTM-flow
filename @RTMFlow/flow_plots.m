function obj = flow_plots(obj)

for i = 1:length(obj.pressure_gradients)
    figure(1);
    set(gcf,'position',[300,300,1500,1000])
    subplot(2,2,1)
    pdeplot(obj.Delaunay_mesh_class.nodes',obj.Delaunay_mesh_class.elements', ...
        XYData=log(obj.physics_class.permeability),XYStyle='flat',ColorMap='jet',Mesh='on')
    axis equal
    xlim([0,1]); ylim([0,1]);
    set(gca, 'YDir', 'normal');
    colormap jet;
    colorbar;
    xlabel('x'); ylabel('y');
    title("Log-permeability")
    
    subplot(2,2,2)
    pdeplot(obj.Delaunay_mesh_class.nodes',obj.Delaunay_mesh_class.elements', ...
        XYData=obj.pressures(:,i),XYStyle='interp',ColorMap='jet',Mesh='on')
    axis equal
    xlim([0,1]); ylim([0,1]);
    set(gca, 'YDir', 'normal');
    colormap jet;
    colorbar;
    xlabel('x'); ylabel('y');
    title("Pressure")

    subplot(2,2,3)
    pdeplot(obj.Delaunay_mesh_class.nodes',obj.Delaunay_mesh_class.elements', ...
        XYData=obj.filling_factors(:,i),XYStyle='interp',ColorMap='jet',Mesh='on')
    axis equal
    xlim([0,1]); ylim([0,1]);
    set(gca, 'YDir', 'normal');
    colormap jet;
    colorbar;
    xlabel('x'); ylabel('y');
    title("Filling factors")

    subplot(2,2,4)
    pdeplot(obj.Delaunay_mesh_class.nodes',obj.Delaunay_mesh_class.elements', ...
        'flowdata',-obj.pressure_gradients{i})
    axis equal
    xlim([0,1]); ylim([0,1]);
    set(gca, 'YDir', 'normal');
    colormap jet
    colorbar;
    xlabel('x'); ylabel('y');
    title("Negative pressure gradient")

    sgtitle("Time elapsed: " + num2str(obj.times(i)) + "s.")
end

end