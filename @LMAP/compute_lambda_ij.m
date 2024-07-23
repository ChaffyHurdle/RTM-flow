function lambda_ij = compute_lambda_ij(obj,i,j,stiffness_matrices,active_matrix,dirichlet_matrix)

% Various init variables
nodes = obj.mesh_class.nodes;
lambda_ij = zeros(length(nodes),length(obj.RTMflow_class.times));
inlet_nodes = obj.RTMflow_class.pressure_class.is_inlet;
t_j = obj.RTMflow_class.observation_times(j);
x_i = obj.RTMflow_class.sensor_locs_on_mesh(i,:);
x_i_ind = obj.RTMflow_class.sensor_inds_on_mesh(i);

% Solve adjoint backwards in time!
for t = length(obj.RTMflow_class.pressure_gradients):-1:1
    
    % Various data from the forward solve
    stiffness_matrix = stiffness_matrices{t};
    active_nodes = active_matrix(:,t); % nodes in D(t)
    dirichlet_nodes = dirichlet_matrix(:,t); % inlet, front and beyond front
    free = active_nodes & ~dirichlet_nodes; % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free
    is_moving_boundary = dirichlet_nodes ...
                   & active_nodes ...
                   & ~inlet_nodes; % nodes that are active, Dirichlet but not inlet
    
    % Set Dirichlet boundary conditions
    lambda = zeros(length(nodes),1);
    lambda(inlet_nodes) = 0.0;
    lambda(is_moving_boundary) = obj.f(obj.RTMflow_class.times(t))*nodes(is_moving_boundary,2)*10000;
    
    load_vector = zeros(length(nodes),1);
    is_sensor_i_active = active_nodes(x_i_ind);
    load_vector(x_i_ind) = obj.delta_t(t_j, obj.RTMflow_class.times(t))*is_sensor_i_active;

%     load_vector = obj.delta_t(t_j, obj.RTMflow_class.times(t))*...
%         obj.delta_x(x_i, nodes);
    b_free = load_vector(free) - ...
                      stiffness_matrix(free,fixed)*lambda(fixed);
    A_free = stiffness_matrix(free,free);

    %% solve matrix equation for pressure
    lambda(free) = A_free\b_free;
    lambda_ij(:,t) = lambda;
end
end