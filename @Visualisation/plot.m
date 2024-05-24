function plot(obj,cvfem_class) 

%% plots moving volume
if obj.is_plotting_volume
    figure(1)
    clf;
    fluid_elements = cvfem_class.active_elements;
    fluid_elements(fluid_elements>0) = 1;

    pdeplot(cvfem_class.mesh_class.nodes',...
            cvfem_class.mesh_class.elements', ...
            XYData=fluid_elements,ColorMap="jet",Mesh="on")
    axis equal
    caxis([0 1]);
    title(["plot of fluid position at time = " ...
                                    num2str(cvfem_class.time)])
end

%% plots pressure solution 
if obj.is_plotting_pressure
    figure(2)
    clf
    
    pressure = cvfem_class.pressure_class.pressure;

    axis equal
    axis auto
    pdeplot(cvfem_class.mesh_class.nodes',...
            cvfem_class.mesh_class.elements', ...
            XYData=pressure,ZData = pressure,...
            ColorMap="jet")
    title(["plot of pressure at time = " ...
                                    num2str(cvfem_class.time)])
end

end