function obj = plotter(obj,u,C,h)

% Plots ongoing LMAP iterations

figure(1)

% Truth
subplot(2,2,1)
c_min = min(min(obj.inverse_class.u_true),min(u));
c_max = max(max(obj.inverse_class.u_true),max(u));
pdeplot(obj.inverse_class.fwd_mesh.nodes',...
        obj.inverse_class.fwd_mesh.elements', ...
        XYData=obj.inverse_class.u_true, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([c_min,c_max])
title("$u^\dagger$",'Interpreter','latex')
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off

% u_k
subplot(2,2,2)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',XYData=u, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([c_min,c_max])
title("$u_{map}$",'interpreter','latex')
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off

% h (step taken)
subplot(2,2,3)
pdeplot(obj.mesh_class.nodes',...
        obj.mesh_class.elements', ...
        XYData=h, ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off

% Variance
title("$h$",'Interpreter','latex')
subplot(2,2,4)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',XYData=diag(C), ...
        XYStyle='interp',ColorMap="jet",Mesh="off")
clim([0,obj.inverse_class.matern_var])
hold on
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
hold off
title("$C_{map}$",'interpreter','latex')


drawnow

end
