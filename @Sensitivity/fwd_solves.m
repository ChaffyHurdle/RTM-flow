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

end