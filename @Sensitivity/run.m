function obj = run(obj)

% Mesh
elements = obj.mesh_class.elements;
nodes = obj.mesh_class.nodes;
element_areas = obj.mesh_class.element_areas;
shape_fun_grads = obj.RTMflow_u.pressure_class.shape_fun_gradients;
num_nodes = length(nodes);
ntimes = length(obj.RTMflow_u.times);

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
%     if t > 50 && t < 60
%         figure(1)
%         subplot(1,2,1)
%         pdeplot(nodes',elements',XYData = double(obj.RTMflow_u.all_active_elements(:,t)),XYStyle='flat',ColorMap="jet",Mesh="on")
%         hold on
%         plot(nodes(is_moving_boundary,1),nodes(is_moving_boundary,2),'wo')
%         hold off
%         subplot(1,2,2)
%         pdeplot(nodes',elements',XYData = pressures(:,t),XYStyle='flat',ColorMap="jet",Mesh="on")
%         hold on
%         plot(nodes(is_moving_boundary,1),nodes(is_moving_boundary,2),'wo')
%         hold off
%         input('')
%     end

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
    p_tilde_t(is_moving_boundary) = max(p_tilde_t(is_moving_boundary),0);

%     figure(1)
%     subplot(1,2,1)
%     pdeplot(nodes',...
%             elements', ...
%             XYData=p_tilde_t,XYStyle='flat',ColorMap="jet",Mesh="on")
%     hold on
%     scatter(nodes(is_moving_boundary,1),nodes(is_moving_boundary,2),'w.')
%     hold off
%     subplot(1,2,2)
%     pdeplot(nodes',...
%             elements', ...
%             XYData=obj.bndry_cond.*is_moving_boundary,XYStyle='flat',ColorMap="jet",Mesh="on")
%     hold on
%     scatter(nodes(is_moving_boundary,1),nodes(is_moving_boundary,2),'w.')
%     hold off

    % Save and update boundary condition
    p_tilde(:,t) = p_tilde_t;

end

obj.p_tilde = p_tilde;

p_tilde_ob_times = [];
is_moving_boundary_ob_times = [];
is_moving_boundary_ob_times_u_plus_h = [];
times_u = obj.RTMflow_u.times;
times_u_plus_h = obj.RTMflow_u_plus_h.times;

for i = 1:length(obj.darcy_class_u.observation_times)
    ob_time = obj.darcy_class_u.observation_times(i);
    idx_u_after = find(times_u > ob_time, 1);
    idx_u_plus_h_after = find(times_u_plus_h > ob_time, 1);
    dt_u = times_u(idx_u_after)-times_u(idx_u_after-1);

    interpolated_p_tilde = p_tilde(:,idx_u_after-1) + ...
                ((ob_time - times_u(idx_u_after-1))/dt_u)*(p_tilde(:,idx_u_after) - p_tilde(:,idx_u_after-1));
    is_moving_boundary_ob_times = [is_moving_boundary_ob_times obj.is_moving_boundary(:,idx_u_after-1)];
    is_moving_boundary_ob_times_u_plus_h = [is_moving_boundary_ob_times_u_plus_h obj.is_moving_boundary_u_plus_h(:,idx_u_plus_h_after-1)];
    p_tilde_ob_times = [p_tilde_ob_times interpolated_p_tilde];
end

obj.p_tildes = p_tilde_ob_times;
obj.is_moving_boundary_ob_times = is_moving_boundary_ob_times;
obj.is_moving_boundary_ob_times_u_plus_h = is_moving_boundary_ob_times_u_plus_h;

num_ob_times = length(obj.darcy_class_u.observation_times);

figure
for i = 1:num_ob_times
    max_diff = max(max(abs(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*obj.active_nodes_u(:,i)),...
        max(abs(obj.pressures_u_plus_h(:,i) - obj.pressures_u(:,i) - obj.p_tildes(:,i)).*obj.active_nodes_u(:,i)));

    % |p_{u+h} - p_{u}|
    subplot(2,num_ob_times,i)
    pdeplot(nodes',elements', ...
        XYData=abs(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([0,max_diff])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$|p_{u+h} - p_u|$", 'interpreter', 'latex')

    % |p_{u+h} - (p_u + \tilde{p})|
    subplot(2,num_ob_times,i+num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=abs(obj.pressures_u_plus_h(:,i) - obj.pressures_u(:,i) - obj.p_tildes(:,i)), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([0,max_diff])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$|p_{u+h} - (p_u + \tilde{p})|$", 'interpreter', 'latex')
end

figure
for i = 1:num_ob_times
    subplot(4,num_ob_times,i)
    max_cbar = max(max(obj.p_tildes(:,i)),max(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)));
    min_cbar = min(min(obj.p_tildes(:,i)),min(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)));

    pdeplot(nodes',elements', ...
        XYData=(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*obj.active_nodes_u(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$p_{u+h}-p_u$",'interpreter','latex')

    subplot(4,num_ob_times,i+num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=obj.p_tildes(:,i).*obj.active_nodes_u(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$\tilde{p}$",'interpreter','latex')

    subplot(4,num_ob_times,i+2*num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*is_moving_boundary_ob_times(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$p_{u+h}-p_u$ on boundary",'interpreter','latex')

    subplot(4,num_ob_times,i+3*num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=obj.p_tildes(:,i).*is_moving_boundary_ob_times(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
%     hold on
%     plot(nodes(boolean(is_moving_boundary_ob_times(:,i)),1),nodes(boolean(is_moving_boundary_ob_times(:,i)),2),'w.')
%     hold off
    title("$\tilde{p}$ on boundary",'interpreter','latex')
end

figure
for i = 1:num_ob_times
    plot(i, sum( (obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i) - obj.p_tildes(:,i)).^2 )/sum(obj.active_nodes_u(:,i)),'k*' )
    hold on
end
hold off

obj.p_tildes_at_sensors = zeros(size(obj.RTMflow_u.pressure_data));
for i = 1:width(obj.p_tildes)
    p_tilde = obj.p_tildes(:,i);
    p_tilde_at_sensors = interpolate_pressures(p_tilde,obj.RTMflow_u,obj.mesh_class);
    obj.p_tildes_at_sensors(:,i) = p_tilde_at_sensors;
end

end

%% The following functions are designed to interpolate pressure at sensor locations
function pressure_at_sensors = interpolate_pressures(p_tilde,obj,mesh)

    sensor_element_inds = obj.sensor_element_inds;
    sensor_locs = obj.physics_class.sensor_locs;
    pressure_at_sensors = zeros(1,length(sensor_locs));

    for i = 1:length(sensor_locs)
        sensor_loc = sensor_locs(i,:);
        nodes_surrounding_inds = mesh.elements(sensor_element_inds(i),:);
        nodes_surrounding = mesh.nodes(nodes_surrounding_inds,:);
        [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, sensor_loc);
        P1 = p_tilde(nodes_surrounding_inds(1));
        P2 = p_tilde(nodes_surrounding_inds(2));
        P3 = p_tilde(nodes_surrounding_inds(3));
        pressure_at_sensors(i) = lambda1 * P1 + lambda2 * P2 + lambda3 * P3;
    end
end

% Auxiliary function (previously defined)
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p_tilde)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p_tilde(1) - v3(1)) + (v3(1) - v2(1))*(p_tilde(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p_tilde(1) - v3(1)) + (v1(1) - v3(1))*(p_tilde(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end