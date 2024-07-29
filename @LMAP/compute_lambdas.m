function obj = compute_lambdas(obj)
    
    %% General comments
    disp("Solving adjoint, representer and linearised state equations...")
    % This code is written to optimise data transfer to all workers.
    % Most notably, we never want to pass obj into the worker function since
    % it contains a lot of (mostly unnecessary) data. This can make this
    % function rather messy.

    % Each worker carries out the following routine: 
    % i -> (lambda_i,kappa_i) -> R_i -> (\tilde{p}_i,\tilde{\Upsilon}_i).
    % This is all done in one pass to avoid opening and closing parfor
    % loops, which would require sending data to workers once more.

    %% Various init variables

    nobs = obj.physics_class.nobservations;
    nsensors = obj.physics_class.nsensors;
    
    % Node data
    nodes = obj.mesh_class.nodes;
    num_nodes = obj.mesh_class.num_nodes;
    inlet_nodes = obj.RTMflow_class.pressure_class.is_inlet;
    sensor_locs = obj.RTMflow_class.sensor_locs_on_mesh;
    sensor_locs_inds = obj.RTMflow_class.sensor_inds_on_mesh;
    
    % Element data
    elements = obj.mesh_class.elements;
    num_elements = obj.mesh_class.num_elements;
    

    times = obj.RTMflow_class.times;
    diff_times = diff(times);
    num_times = length(times);
    ob_times = obj.physics_class.observation_times;
    ob_times_inds = zeros(1,length(ob_times));
    for i = 1:length(ob_times)
        ob_times_inds(i) = find(obj.RTMflow_class.times > ob_times(i),1)-1;
    end

    [i_vec, j_vec] = meshgrid(1:nsensors, 1:nobs);
    i_vec = reshape(i_vec, [], 1);
    j_vec = reshape(j_vec, [], 1);

    lambda_mat = zeros(num_nodes, num_times, nobs * nsensors);
    lambda_grad_mat = zeros(num_elements,2,num_times,nobs * nsensors);
    representers = zeros(num_elements,nsensors*nobs);

    % Broadcast large objects to share between pool
    disp("Moving stiffness to pool")
    broadcast_stiffness_matrices = parallel.pool.Constant(obj.RTMflow_class.stiffness_matrices);
    disp("Moving active nodes to pool")
    active_matrix = parallel.pool.Constant(obj.RTMflow_class.active_nodes); % nodes in D(t)
    disp("Moving dirichlet nodes to pool")
    dirichlet_matrix = parallel.pool.Constant(obj.RTMflow_class.Dirichlet_nodes); % inlet, front and beyond front
    disp("Moving shape func gradients to pool")
    shape_func_grads = parallel.pool.Constant(obj.pressure_class.shape_fun_gradients);
    deltat = obj.delta_t;
    disp("Moving grad pressures to pool")
    grad_pressures = parallel.pool.Constant(obj.RTMflow_class.pressure_gradients);
    permeability = obj.physics_class.permeability;
    disp("Moving active_elements to pool")
    active_elements = parallel.pool.Constant(obj.RTMflow_class.all_active_elements);

    % Compute each adjoint for pressure (along with its gradient)
    parfor k = 1:nobs * nsensors

        i = i_vec(k);
        j = j_vec(k);
        disp([i,j])

        % Retrieve shared data
        stiffness_matrix_ij = broadcast_stiffness_matrices.Value;
        active_matrix_ij = active_matrix.Value;
        dirichlet_matrix_ij = dirichlet_matrix.Value;
        shape_fun_grads_ij = shape_func_grads.Value;
        
        % Access other required variables
        t_j = ob_times(j);
        x_i_ind = sensor_locs_inds(i);
        ob_time_ind = ob_times_inds(j);
        
        % Compute lambda_ij
        lambda_ij = compute_lambda_ij(x_i_ind,t_j,ob_time_ind, ...
            nodes,inlet_nodes,times,deltat, ...
            stiffness_matrix_ij,active_matrix_ij,dirichlet_matrix_ij);
        
        % Store result in lambda_mat
        lambda_mat(:, :, k) = lambda_ij;

        % Compute grad lambda_ij
        grad_lambda_ij = compute_grad_lambda_ij(lambda_ij,shape_fun_grads_ij, ...
    elements,num_elements,num_times);
        lambda_grad_mat(:,:,:,k) = grad_lambda_ij;
        
        % Compute representer_ij
        grad_pressures_i = grad_pressures.Value;
        active_elements_i = active_elements.Value;
        R_i = compute_representer_i(times,diff_times, ...
            grad_lambda_ij,grad_pressures_i,active_elements_i,permeability);
        representers(:,k) = R_i;

    end

    % Set data
    obj.lambdas = lambda_mat;
    obj.grad_lambdas = lambda_grad_mat;
    obj.representers = representers;
end



%% Worker function for lambda_ij
function lambda_ij = compute_lambda_ij(x_i_ind,t_j,ob_time_ind, ...
    nodes,inlet_nodes,times,deltat, ...
    stiffness_matrices,active_matrix,dirichlet_matrix)

% Init
lambda_ij = zeros(length(nodes),length(times));
%del_t = 3*max(diff(times));
del_t = 10e-11;

% Solve adjoint backwards in time!
% if active_matrix(x_i_ind,ob_time_ind)
for t = length(times):-1:1
    
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
    lambda(is_moving_boundary) = 0.0;
    
    load_vector = zeros(length(nodes),1);
    is_sensor_i_active = active_nodes(x_i_ind);
    load_vector(x_i_ind) = deltat(t_j,times(t),del_t)*is_sensor_i_active;

%     load_vector = obj.delta_t(t_j, obj.RTMflow_class.times(t))*...
%         obj.delta_x(x_i, nodes);
    b_free = load_vector(free) - ...
                      stiffness_matrix(free,fixed)*lambda(fixed);
    A_free = stiffness_matrix(free,free);

    %% solve matrix equation for pressure
    lambda(free) = A_free\b_free;
    lambda_ij(:,t) = lambda;
end
%end

end

%% Worker function for grad lambda_ij
function grad_lambda_ij = compute_grad_lambda_ij(lambda_ij,shape_fun_grads, ...
    elements,num_elements,num_times)

grad_lambda_ij = zeros(num_elements,2,num_times);

for i = 1:num_times

    lambda_ij_time = lambda_ij(:,i);

    for j = 1:num_elements 
        
        local_lambda = lambda_ij_time(elements(j,:));
        grad_phi = shape_fun_grads{j};
    
        grad_lambda_ij(j,:,i) = (grad_phi * local_lambda)';
    end
end

end


%% Worker function for representer_i
function R_i = compute_representer_i(times,diff_times, ...
        grad_lambda_i,grad_pressures_i,active_elements_i,permeability)
    R_i = zeros(1,length(grad_lambda_i));

    for i = 1:(length(times)-1)
        f_i = active_elements_i(:,i) .* sum(grad_lambda_i(:,:,i) .* grad_pressures_i{i},2);
        f_iplus1 = active_elements_i(:,i+1) .* sum(grad_lambda_i(:,:,i+1) .* grad_pressures_i{i+1},2);
        R_i = R_i + diff_times(i) * (f_i + f_iplus1)';
    end
    R_i = - R_i .* exp(permeability);

    R_i(abs(R_i)<1e-10)=0;
end