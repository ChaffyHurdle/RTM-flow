function [bndry_cond,v_h_new] = compute_bndry_condition(t,v_h_old,grad_ptilde_ij,h,obj_data)


% Various reused objects
mesh_class = obj_data.mesh_class;
dt_vec = diff(obj_data.RTMflow_class.times);
viscosity = obj_data.RTMflow_class.physics_class.viscosity;
porosity = obj_data.RTMflow_class.physics_class.porosity;
pressure_gradient_t = obj_data.RTMflow_class.pressure_gradients{t};
pressure_gradient_t_minus_1 = obj_data.RTMflow_class.pressure_gradients{t-1};
is_moving_boundary_t = obj_data.RTMflow_class.moving_boundary(:,t);
is_moving_boundary_tminus1 = obj_data.RTMflow_class.moving_boundary(:,t-1);
v_h_new = v_h_old;
u = obj_data.u;
bndry_cond = zeros(1,mesh_class.num_nodes);

%% Edges
if t > 4
    moving_boundary_inds = find(is_moving_boundary_t);

    edge_data = obj_data.RTMflow_class.edge_data{t};

    % Find new active elements, assign value based on nearest element at t-1
    new_active_elements = find(obj_data.RTMflow_class.all_new_active_elements(:,t));
    prev_nodes = is_moving_boundary_tminus1 | obj_data.RTMflow_class.pressure_class.is_inlet;
    prev_nodes_inds = find(prev_nodes);

    prev_elems = zeros(mesh_class.num_elements,1);
    for i = 1 : length(prev_nodes_inds)
        candidate = prev_nodes_inds(i);
        for j = 2 : obj_data.RTMflow_class.Voronoi_mesh_class.has_node_i(candidate,1)+1
            elem = obj_data.RTMflow_class.Voronoi_mesh_class.has_node_i(candidate,j);
            if obj_data.RTMflow_class.all_active_elements(elem,t-1) >= 0.5
                prev_elems(elem) = 1;
            end
            %candidate_elem(elem) = 1;
        end
    end

    prev_elems = find(prev_elems);
    nearest_inds = zeros(mesh_class.num_elements,1);
    for i = 1:length(new_active_elements)
        elem = new_active_elements(i);
        centroid = mesh_class.centroids(elem,:);
        active_centroids = mesh_class.centroids(prev_elems,:);
        distances = sqrt(sum((active_centroids - centroid).^2, 2));
        [~, idx] = min(distances);
        nearest_ind = prev_elems(idx);
        nearest_inds(elem) = nearest_ind;
        %obj.v_h(elem,:) = obj.v_h(nearest_ind,:);

        [~,b] = sort(distances);
        v_h_new(elem,:) = mean(v_h_old(prev_elems(b(1:3)),:));
    end
    
    
    % Triangles on boundary
    candidate_elem = zeros(mesh_class.num_elements,1);
    candidate_elem(edge_data(:,1)) = 1;
    candidate_elem_inds = find(candidate_elem);

    bnd_vals = zeros(mesh_class.num_elements,1);
    
    for i = 1:length(candidate_elem_inds)
    
        elem_index = candidate_elem_inds(i);
    
        %% extract linearised velocity in element centre
        if ismember(elem_index, new_active_elements)
            nearest_element_ind = nearest_inds(elem);
            dvh_dt = - (1/(viscosity*porosity)) * ( exp(u(elem_index)) * h(elem_index) * pressure_gradient_t_minus_1(nearest_element_ind,:)' ...
                              + exp(u(elem_index)) * grad_ptilde_ij(nearest_element_ind,:)');
        else
            dvh_dt = - (1/(viscosity*porosity)) * ( exp(u(elem_index)) * h(elem_index) * pressure_gradient_t_minus_1(elem_index,:)' ...
                              + exp(u(elem_index)) * grad_ptilde_ij(elem_index,:)' );
        end
        
        v_h_new(elem_index,:) = v_h_new(elem_index,:) + dt_vec(t-1)*dvh_dt';
    
        %% Compute transfer into each CV
        local_linearised_flow_rate = local_flux_tri(edge_data(i,4:5),v_h_new(elem_index,:)');
        local_pressure_grad = local_flux_tri(edge_data(i,4:5),pressure_gradient_t(elem_index,:)');
    
        %% Need to track v_h as it depends on time
        bnd_vals(elem_index) = - local_linearised_flow_rate*local_pressure_grad;
    
    end
    %obj.v_h(obj.RTMflow_u.pressure_class.is_inlet,:) = 0;

    for i = 1:length(moving_boundary_inds)
        node_ind = moving_boundary_inds(i);
        rows_involving_node = any(edge_data(:,2:3) == node_ind,2);
        bndry_cond(node_ind) = mean(bnd_vals(edge_data(rows_involving_node,1)));
    end
    
end

end

function qn = local_flux_tri(normal_vec,v)

qn = normal_vec*v;

end
