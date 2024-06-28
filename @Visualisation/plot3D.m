function plot3D(obj,RTM_class) 

%% plots moving volume
if obj.is_plotting_volume
    figure(1)
    clf;
    volume = zeros(RTM_class.Delaunay_mesh_class.num_nodes,1);
    volume(RTM_class.pressure_class.is_Dirichlet) = 1;

    pdeplot3D(RTM_class.Delaunay_mesh_class.nodes',...
            RTM_class.Delaunay_mesh_class.elements', ...
            ColorMapData=volume)

    axis equal
    clim([0 1]);
    title(["plot of fluid position at time = " ...
                                    num2str(RTM_class.time)])

    if obj.is_animate_volume 
        exportgraphics(gca,"volume.gif","Append",true)
    end

end

if obj.is_plotting_pressure
figure(2)
clf;

pdeplot3D(RTM_class.Delaunay_mesh_class.nodes',...
        RTM_class.Delaunay_mesh_class.elements', ...
        ColorMapData=RTM_class.pressure_class.pressure,mesh="on")

axis equal
title(["plot of fluid position at time = " ...
                                num2str(RTM_class.time)])

if obj.is_animate_volume 
    exportgraphics(gca,"volume.gif","Append",true)
end
end


end