function obj = run(obj)

% Various data/inits
elements = obj.mesh_class.elements;
nodes = obj.mesh_class.nodes;
element_areas = obj.mesh_class.element_areas;
shape_fun_grads = obj.RTMflow_u.pressure_class.shape_fun_gradients;
num_nodes = length(nodes);
ntimes = length(obj.RTMflow_u.times);

p_tilde = zeros(num_nodes,ntimes);
stiffness_left = obj.RTMflow_u.stiffness_matrices;
stiffness_right = sparse(num_nodes,num_nodes,10*num_nodes);
active_matrix = obj.RTMflow_u.active_nodes;
all_new_active_elements = [zeros(length(elements),1) obj.RTMflow_u.all_new_active_elements];
dirichlet_matrix = obj.RTMflow_u.Dirichlet_nodes;
inlet_nodes = obj.RTMflow_u.pressure_class.is_inlet;
pressures = obj.RTMflow_u.pressures;

for t = 2:ntimes

    stiffness_left_t = stiffness_left{t-1};
    active_nodes = boolean(active_matrix(:,t-1)); % nodes in D(t)
    dirichlet_nodes = boolean(dirichlet_matrix(:,t-1)); % inlet, front and beyond front
    free = active_nodes & ~dirichlet_nodes; % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free
    is_moving_boundary = dirichlet_nodes ...
                   & active_nodes ...
                   & ~inlet_nodes; % nodes that are active, Dirichlet but not inlet
    
    % Set Dirichlet boundary conditions
    p_tilde_t = zeros(num_nodes,1);
    p_tilde_t(inlet_nodes) = 0.0;
    p_tilde_t(is_moving_boundary) = 0.0;

    % Set RHS
    new_elements = find(all_new_active_elements(:,t));
    for i = 1 : length(new_elements)
    
        %% element geometry
        tri_nodes = elements(new_elements(i),:);
        area = element_areas(new_elements(i));
        
        %% local FEM stiffness matrix
        grad_phi = shape_fun_grads{new_elements(i)};
        A_local = exp(obj.u(new_elements(i)))*obj.h(new_elements(i)) *(grad_phi'*grad_phi)*area;
        
        %% local to global mapping
        stiffness_right(tri_nodes,tri_nodes) = stiffness_right(tri_nodes,tri_nodes) + A_local;
    end
    b_free = - stiffness_right(free,active_nodes)*pressures(active_nodes,t) - stiffness_left_t(free,fixed)*p_tilde_t(fixed);
    A_free = stiffness_left_t(free,free);
    p_tilde_t(free) = A_free\b_free;
    p_tilde(:,t) = p_tilde_t;
end

obj.p_tilde = p_tilde;

end