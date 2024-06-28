function plot3D(obj,RTM_class) 

%% plots moving volume
if obj.is_plotting_volume
    figure(1)
    clf;
    fluid_elements = RTM_class.active_elements;
    fluid_elements(fluid_elements>0) = 1;

    pdeplot3D(RTM_class.Delaunay_mesh_class.nodes',...
            RTM_class.Delaunay_mesh_class.elements', ...
            ColorMapData=fluid_elements)

    axis equal
    clim([0 1]);
    title(["plot of fluid position at time = " ...
                                    num2str(RTM_class.time)])

    if obj.is_animate_volume 
        exportgraphics(gca,"volume.gif","Append",true)
    end

end