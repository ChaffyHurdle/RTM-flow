function obj = update_boundary_condition(obj,t)

% Various variable sets for brevity
mesh_class = obj.mesh_class;
dt_vec = diff(obj.RTMflow_u.times);
viscosity = obj.RTMflow_u.physics_class.viscosity;
porosity = obj.RTMflow_u.physics_class.porosity;
normal_vec = obj.RTMflow_u.Voronoi_mesh_class.volume_outflow_vectors;
pressure_gradient_t = obj.RTMflow_u.pressure_gradients{t};

% Find nodes associated with inlet, and the boundary at [t,t+1]
is_moving_boundary = obj.is_moving_boundary(:,t) | obj.is_moving_boundary(:,t+1);
inlet_nodes = obj.RTMflow_u.pressure_class.is_inlet;
dirichlet_nodes = inlet_nodes | is_moving_boundary;
obj.dirichlet_nodes_matrix = [obj.dirichlet_nodes_matrix dirichlet_nodes];
dirichlet_inds = find(dirichlet_nodes);

% Find all elements connected to Dirichlet nodes
candidate_elem = zeros(length(mesh_class.elements),1);
for i = 1 : length(dirichlet_inds)
    candidate = dirichlet_inds(i);
    for j = 2 : obj.RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,1)+1
        elem = obj.RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,j);
        if obj.RTMflow_u.all_active_elements(elem,t) >= 0.5 || obj.RTMflow_u.all_active_elements(elem,t+1) >= 0.5
            candidate_elem(elem) = 1;
        end
        %candidate_elem(elem) = 1;
    end
end

candidate_elem = find(candidate_elem);

% Set pressure gradient to 0
pressure_grad = zeros(length(mesh_class.nodes),1);
dvh_dt = zeros(mesh_class.num_nodes,1);

% For each element, compute the transfer of v_h and \nabla p between CVs
for i = 1:length(candidate_elem)

    elem_index = candidate_elem(i);

    %% extract local element properties
    element = mesh_class.elements(elem_index,:);

    %% extract linearised velocity in element centre
    tilde_vi = - ( exp(obj.u(elem_index)) * obj.h(elem_index) * pressure_gradient_t(elem_index,:)' ...
                      + exp(obj.u(elem_index)) * obj.grad_p_tilde(elem_index,:)' )/(viscosity*porosity);

    %% Compute transfer into each CV
    local_linearised_flow_rate = local_flux_tri(normal_vec(elem_index,:),tilde_vi,obj.RTMflow_u.physics_class);
    local_pressure_grad = local_flux_tri(normal_vec(elem_index,:),pressure_gradient_t(elem_index,:)',obj.RTMflow_u.physics_class);

    %% Need to track v_h as it depends on time
    dvh_dt(element) = dvh_dt(element) + local_linearised_flow_rate;
    pressure_grad(element) = pressure_grad(element) + local_pressure_grad;

end
obj.v_h = obj.v_h + dvh_dt*dt_vec(t)./obj.RTMflow_u.Voronoi_mesh_class.volume_measures;
obj.bndry_cond = - obj.v_h .* pressure_grad;

%% function to compute outflow of standard elements
function qn = local_flux_tri(normal_vec,v,darcy_class)
qn = zeros(3,1);

qn(1) = normal_vec(:,[1 2])*v;
qn(2) = normal_vec(:,[3 4])*v;
qn(3) = normal_vec(:,[5 6])*v;

qn = qn*darcy_class.thickness;

end

end