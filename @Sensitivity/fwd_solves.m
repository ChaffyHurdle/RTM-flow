function obj = fwd_solves(obj)

mesh = obj.mesh_class;

u = obj.u;
h = obj.h;
darcy_class_u = obj.darcy_class_u;
my_pressure_u = Pressure(mesh,darcy_class_u);
RTMflow_u = RTMFlow(mesh,darcy_class_u,my_pressure_u);
RTMflow_u = RTMflow_u.run();
obj.RTMflow_u = RTMflow_u;

n = 5;
RTMflow_u_plus_h_saves = cell(1,n);

parfor i = 1:n
    disp(i)
    u_plus_h = u + 2^(-i)*h;
    darcy_class_u_plus_h = darcy_class_u;
    darcy_class_u_plus_h.permeability = exp(u_plus_h);
    my_pressure_u_plus_h = Pressure(mesh,darcy_class_u_plus_h);
    RTMflow_u_plus_h = RTMFlow(mesh,darcy_class_u_plus_h,my_pressure_u_plus_h);
    RTMflow_u_plus_h = RTMflow_u_plus_h.run();
    RTMflow_u_plus_h_saves{i} = RTMflow_u_plus_h; 
end

obj.RTMflow_u_plus_h = RTMflow_u_plus_h_saves;

end