function obj = run(obj)

% Runs EKI.run_t(t) method for each time.

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

% Save data
obj.ueki_seq = ueki_seq;
obj.Ceki_seq = Ceki_seq;
obj.timer_seq = timer_seq;