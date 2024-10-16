function p_tilde_ij = compute_ptilde_ij(R_ij,obj_data)

% Mesh
elements = obj_data.mesh_class.elements;
element_areas = obj_data.mesh_class.element_areas;
shape_fun_grads = obj_data.RTMflow_class.pressure_class.shape_fun_gradients;
num_nodes = obj_data.mesh_class.num_nodes;
num_elements = obj_data.mesh_class.num_elements;
ntimes = length(obj_data.RTMflow_class.times);
u = obj_data.u;
h = R_ij/(1+obj_data.alpha);

% Forward simulation data
stiffness_left = obj_data.RTMflow_class.stiffness_matrices;
stiffness_right = spalloc(num_nodes,num_nodes,10*num_nodes);
active_matrix = obj_data.RTMflow_class.active_nodes;
dirichlet_matrix = obj_data.RTMflow_class.Dirichlet_nodes;
inlet_nodes = obj_data.RTMflow_class.pressure_class.is_inlet;
pressures = obj_data.RTMflow_class.pressures;

% Initialise
v_h = zeros(num_elements,2);
p_tilde_ij = zeros(num_nodes,ntimes);
grad_ptilde_ij = zeros(num_elements,2);

% Loop over forward simulation
for t = 2:ntimes-1
    % Left side of the variational form is the same as the forward problem
    stiffness_left_t = stiffness_left{t};
    
    % Various collections of nodes
    active_nodes = logical(active_matrix(:,t)); % nodes in D(t) (inc. edges)
    dirichlet_nodes = logical(dirichlet_matrix(:,t)); % inlet, front and beyond front
    free = active_nodes & ~dirichlet_nodes; % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free (i.e. Dirichlet)
    is_moving_boundary = obj_data.RTMflow_class.moving_boundary(:,t); % nodes that are active, Dirichlet but not inlet

    % Set Dirichlet boundary conditions
    p_tilde_t = zeros(num_nodes,1);
    [bndry_cond, v_h] = compute_bndry_condition(t,v_h,grad_ptilde_ij,h,obj_data);
    p_tilde_t(is_moving_boundary) = bndry_cond(is_moving_boundary);
    p_tilde_t(inlet_nodes) = 0.0;
    
    % Find RHS stiffness in variational form
    new_elements = find(obj_data.RTMflow_class.all_new_active_elements(:,t));
    for i = 1 : length(new_elements)
    
        %% element geometry
        tri_nodes = elements(new_elements(i),:);
        area = element_areas(new_elements(i));
        
        %% local FEM stiffness matrix
        grad_phi = shape_fun_grads(:,:,new_elements(i));
        A_local = (exp(u(new_elements(i)))*h(new_elements(i))*grad_phi)' *grad_phi*area;
        
        %% local to global mapping
        stiffness_right(tri_nodes,tri_nodes) = stiffness_right(tri_nodes,tri_nodes) + A_local;
    end
    % Set up linear system
    b_free = - stiffness_right(free,active_nodes)*pressures(active_nodes,t) - stiffness_left_t(free,fixed)*p_tilde_t(fixed);
    A_free = stiffness_left_t(free,free);
    p_tilde_t(free) = A_free\b_free;
    
    grad_ptilde_ij = compute_grad_ptilde_ij(p_tilde_t,obj_data);
    
    % Save and update boundary condition
    p_tilde_ij(:,t) = p_tilde_t;


end
end


function grad_ptilde_ij = compute_grad_ptilde_ij(ptilde,obj_data)

elements = obj_data.mesh_class.elements;
shape_fun_grads = obj_data.RTMflow_class.pressure_class.shape_fun_gradients;

local_ptildes = ptilde(elements);
grad_ptilde_ij = squeeze(sum(bsxfun(@times, shape_fun_grads, permute(local_ptildes, [3, 2, 1])), 2))';

end