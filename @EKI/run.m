function obj = run(obj)

start_permeability = obj.physics_class.permeability;
start_pressure_class = obj.pressure_class;
start_flow_class = obj.RTMflow_class;

n = obj.physics_class.nobservations;
ueki_seq = zeros(obj.mesh_class.num_elements,n);
Ceki_seq = zeros(obj.mesh_class.num_elements,obj.mesh_class.num_elements,n);
timer_seq = zeros(1,n);

for t = 1:n
    obj = obj.run_t(t);
    
    ueki_seq(:,t) = obj.ueki;
    Ceki_seq(:,:,t) = obj.Ceki;
    timer_seq(t) = obj.timer;

end
obj.ueki_seq = ueki_seq;
obj.Ceki_seq = Ceki_seq;
obj.timer_seq = timer_seq;