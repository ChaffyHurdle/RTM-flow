function obj = fwd_solves(obj)

mesh = obj.mesh_class;

u = obj.u;
h = obj.h;
darcy_class_u = obj.darcy_class_u;
my_pressure_u = Pressure(mesh,darcy_class_u);
RTMflow_u = RTMFlow(mesh,darcy_class_u,my_pressure_u);
RTMflow_u = RTMflow_u.run();
obj.RTMflow_u = RTMflow_u;

u_plus_h = u + h;
darcy_class_u_plus_h = darcy_class_u;
darcy_class_u_plus_h.permeability = exp(u_plus_h);
my_pressure_u_plus_h = Pressure(mesh,darcy_class_u_plus_h);
RTMflow_u_plus_h = RTMFlow(mesh,darcy_class_u_plus_h,my_pressure_u_plus_h);
RTMflow_u_plus_h = RTMflow_u_plus_h.run(); 
obj.RTMflow_u_plus_h = RTMflow_u_plus_h;

obj.is_moving_boundary = RTMflow_u.active_nodes & RTMflow_u.Dirichlet_nodes & ~my_pressure_u.is_inlet;
obj.is_moving_boundary_u_plus_h = RTMflow_u_plus_h.active_nodes & RTMflow_u_plus_h.Dirichlet_nodes & ~my_pressure_u_plus_h.is_inlet;

pressure_u = RTMflow_u.pressures;
times_u = RTMflow_u.times;
pressure_u_plus_h = RTMflow_u_plus_h.pressures;
times_u_plus_h = RTMflow_u_plus_h.times;

pressures_u = [];
pressures_u_plus_h = [];
active_nodes_u = [];
time_inds_u = [];
time_inds_u_plus_h = [];

for i = 1:length(obj.darcy_class_u.observation_times)
    ob_time = obj.darcy_class_u.observation_times(i);
    idx_u = find(times_u > ob_time, 1) - 1;
    time_inds_u =[time_inds_u idx_u];
    idx_u_plus_h = find(times_u_plus_h > ob_time, 1) - 1;
    time_inds_u_plus_h = [time_inds_u_plus_h idx_u_plus_h];
    active_nodes_u = [active_nodes_u RTMflow_u.active_nodes(:,idx_u)];
%     interpolated_p_u = pressure_u(:,idx_u_after-1) + ...
%                 ((ob_time - times_u(idx_u_after-1))/dt_u)*(pressure_u(:,idx_u_after) - pressure_u(:,idx_u_after-1));
%     interpolated_p_u_plus_h = pressure_u_plus_h(:,idx_u_plus_h_after-1) + ...
%                 ((ob_time - times_u_plus_h(idx_u_plus_h_after-1))/dt_u_plus_h)*(pressure_u_plus_h(:,idx_u_plus_h_after) - pressure_u_plus_h(:,idx_u_plus_h_after-1));
    pressures_u = [pressures_u pressure_u(:,idx_u)];
    pressures_u_plus_h = [pressures_u_plus_h pressure_u_plus_h(:,idx_u_plus_h)];
end

obj.pressures_u = pressures_u;
obj.pressures_u_plus_h = pressures_u_plus_h;
obj.time_inds_u = time_inds_u;
obj.time_inds_u_plus_h = time_inds_u_plus_h;
obj.active_nodes_u = active_nodes_u;
obj.RTMflow_u.all_new_active_elements = [zeros(mesh.num_elements,1) obj.RTMflow_u.all_new_active_elements];