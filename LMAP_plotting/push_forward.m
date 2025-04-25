function [pressures,flow_fronts] = push_forward(u_samples, physics_class, ...
    mesh_class)

n = size(u_samples,1);
data_size = physics_class.nobservations * physics_class.nsensors;
pressures = zeros(n,data_size);
flow_fronts = cell(n,5);

parfor i = 1:n
    disp(i)
    K_i = exp(u_samples(i,:)');
    physics_class_i = physics_class;
    physics_class_i.permeability = K_i;
    pressure_class_i = Pressure(mesh_class,physics_class_i);
    RTMflow_class_i = RTMFlow(mesh_class,physics_class_i,pressure_class_i);
    RTMflow_class_i = RTMflow_class_i.run(inf);
    pressures(i,:) = reshape(transpose(RTMflow_class_i.pressure_data),1,[]);

    for j = 1:5
        t_index = find(RTMflow_class_i.times > physics_class.observation_times(j),1)-1;
        flow_fronts{i,j} = RTMflow_class_i.edge_data{t_index};
    end
end
