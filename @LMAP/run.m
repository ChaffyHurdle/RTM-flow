function obj = run(obj)

% Runs LMAP.run_t(t) for t = t_1, ..., t_N

start_permeability = obj.physics_class.permeability;
start_pressure_class = obj.pressure_class;
start_flow_class = obj.RTMflow_class;
start_alpha = obj.alpha;

n = obj.physics_class.nobservations;
umap_seq = zeros(obj.mesh_class.num_elements,n);
Cmap_seq = zeros(obj.mesh_class.num_elements,obj.mesh_class.num_elements,n);
timer_seq = zeros(1,n);

for t = 1:n
    disp(t)
    obj.physics_class.permeability = start_permeability;
    obj.pressure_class = start_pressure_class;
    obj.RTMflow_class = start_flow_class;
    obj.alpha = start_alpha;
    obj.u = obj.inverse_class.u0;

    obj = obj.run_t(t);
    
    umap_seq(:,t) = obj.u_map;
    Cmap_seq(:,:,t) = obj.C_map;
    timer_seq(t) = obj.execution_times(end);

end
obj.umap_seq = umap_seq;
obj.Cmap_seq = Cmap_seq;
obj.timer_seq = timer_seq;