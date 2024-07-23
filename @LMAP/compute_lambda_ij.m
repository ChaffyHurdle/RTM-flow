function lambda_ij = compute_lambda_ij(obj,i,j)

nodes = obj.mesh_class.nodes;
lambda_ij = zeros(length(nodes),length(obj.RTMflow_class.times));
inlet_nodes = obj.RTMflow_class.pressure_class.is_inlet;
t_j = obj.RTMflow_class.observation_times(j);
x_i = obj.RTMflow_class.sensor_locs_on_mesh(i,:);

for t = length(obj.RTMflow_class.pressure_gradients):-1:1

    stiffness_matrix = obj.RTMflow_class.stiffness_matrices{t};
    active_nodes = obj.RTMflow_class.active_nodes(:,t);
    dirichlet_nodes = obj.RTMflow_class.Dirichlet_nodes(:,t);
    free = active_nodes & ~dirichlet_nodes; fixed = ~free;
    is_moving_boundary = dirichlet_nodes ...
                   & active_nodes ...
                   & ~inlet_nodes;

    lambda = zeros(length(nodes),1);
    lambda(inlet_nodes) = 0.0;
    lambda(is_moving_boundary) = obj.f(obj.RTMflow_class.times(t))*nodes(is_moving_boundary,2)*10000;

    load_vector = obj.delta_t(t_j, obj.RTMflow_class.times(t))*...
        obj.delta_x(x_i, nodes);
    b_free = load_vector(free)' - ...
                      stiffness_matrix(free,fixed)*lambda(fixed);
    A_free = stiffness_matrix(free,free);

    %% solve matrix equation for pressure
    lambda(free) = A_free\b_free;
    lambda_ij(:,t) = lambda;
end
end