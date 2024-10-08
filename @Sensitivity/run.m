function obj = run(obj)

% Mesh
elements = obj.mesh_class.elements;
nodes = obj.mesh_class.nodes;
element_areas = obj.mesh_class.element_areas;
shape_fun_grads = obj.RTMflow_u.pressure_class.shape_fun_gradients;
num_nodes = length(nodes);
ntimes = length(obj.RTMflow_u.times);
k = 1;

% Forward simulation data
stiffness_left = obj.RTMflow_u.stiffness_matrices;
stiffness_right = spalloc(num_nodes,num_nodes,10*num_nodes);
active_matrix = obj.RTMflow_u.active_nodes;
all_new_active_elements = obj.RTMflow_u.all_new_active_elements;
dirichlet_matrix = obj.RTMflow_u.Dirichlet_nodes;
inlet_nodes = obj.RTMflow_u.pressure_class.is_inlet;
pressures = obj.RTMflow_u.pressures;

% Initialise
obj.v_h = zeros(obj.mesh_class.num_elements,2);
obj.bndry_cond = zeros(obj.mesh_class.num_nodes,1);
p_tilde = zeros(num_nodes,ntimes);

% Loop over forward simulation
for t = 2:ntimes-1
    % Left side of the variational form is the same as the forward problem
    stiffness_left_t = stiffness_left{t};

    % Various collections of nodes
    active_nodes = boolean(active_matrix(:,t)); % nodes in D(t) (inc. edges)
    dirichlet_nodes = boolean(dirichlet_matrix(:,t)); % inlet, front and beyond front
    free = active_nodes & ~dirichlet_nodes; % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free (i.e. Dirichlet)
    is_moving_boundary = obj.is_moving_boundary(:,t); % nodes that are active, Dirichlet but not inlet

    % Set Dirichlet boundary conditions
    obj.bndry_conds = [obj.bndry_conds obj.bndry_cond];
    p_tilde_t = zeros(num_nodes,1);
    obj = obj.update_boundary_condition3(t);
    p_tilde_t(is_moving_boundary) = obj.bndry_cond(is_moving_boundary);
    p_tilde_t(inlet_nodes) = 0.0;

    % Find RHS stiffness in variational form
    new_elements = find(all_new_active_elements(:,t-1));

    for i = 1 : length(new_elements)
    
        %% element geometry
        tri_nodes = elements(new_elements(i),:);
        area = element_areas(new_elements(i));
        
        %% local FEM stiffness matrix
        grad_phi = shape_fun_grads{new_elements(i)};
        A_local = (exp(obj.u(new_elements(i)))*obj.h(new_elements(i))*grad_phi)' *grad_phi*area;
        
        %% local to global mapping
        stiffness_right(tri_nodes,tri_nodes) = stiffness_right(tri_nodes,tri_nodes) + A_local;
    end
    % Set up linear system
    b_free = - stiffness_right(free,active_nodes)*pressures(active_nodes,t) - stiffness_left_t(free,fixed)*p_tilde_t(fixed);
    A_free = stiffness_left_t(free,free);
    p_tilde_t(free) = A_free\b_free;
    
    obj = obj.compute_grad_p_tilde(p_tilde_t);
 
    if ismember(t,obj.time_inds_u)
        figure(10)
        subplot(1,5,k)
        pdeplot(nodes',...
                elements', ...
                XYData=p_tilde_t.*free,XYStyle='flat',ColorMap="jet",Mesh="off")
        hold on
        %scatter(nodes(is_moving_boundary,1),nodes(is_moving_boundary,2),'w.')
        t_uplush = find(obj.RTMflow_u_plus_h.times > obj.RTMflow_u.times(t),1)-1;
        %scatter(nodes(obj.is_moving_boundary_u_plus_h(:,t_uplush),1),nodes(obj.is_moving_boundary_u_plus_h(:,t_uplush),2),'wo')
        for j = 1:size(obj.all_edge_data{t},1)
            plot([nodes(obj.all_edge_data{t}(j,2),1), nodes(obj.all_edge_data{t}(j,3),1)],...
                [nodes(obj.all_edge_data{t}(j,2),2), nodes(obj.all_edge_data{t}(j,3),2)],'w','LineWidth',2);
            midpoint = (nodes(obj.all_edge_data{t}(j,2),:) + obj.mesh_class.nodes(obj.all_edge_data{t}(j,3),:))/2;
            vdotn = dot(obj.v_h(obj.all_edge_data{t}(j,1),:), obj.all_edge_data{t}(j,4:5));
            plot([midpoint(1), midpoint(1)+vdotn*obj.all_edge_data{t}(j,4)],[midpoint(2),midpoint(2) + vdotn*obj.all_edge_data{t}(j,5)],'w-','LineWidth',1)
        end
        mov_bndry_uplush = obj.is_moving_boundary_u_plus_h(:,t_uplush);
        mov_bdnry_uplush_inds = find(mov_bndry_uplush);
        mov_bdnry_uplush_nodes = nodes(mov_bdnry_uplush_inds,:);
        [~,sorted_inds_increasing] = sort(mov_bdnry_uplush_nodes(:,2));
        for j = 1:length(sorted_inds_increasing)-1
           plot([nodes(mov_bdnry_uplush_inds(sorted_inds_increasing(j)),1), nodes(mov_bdnry_uplush_inds(sorted_inds_increasing(j+1)),1)],...
                [nodes(mov_bdnry_uplush_inds(sorted_inds_increasing(j)),2), nodes(mov_bdnry_uplush_inds(sorted_inds_increasing(j+1)),2)],'m-','LineWidth',2); 
        end
        k = k+1;
        hold off
        xlim([0,1]); ylim([0,1]);
    end

    % Save and update boundary condition
    p_tilde(:,t) = p_tilde_t;

end

obj.p_tilde = p_tilde;

obj = obj.save_sensitivity_data();

end