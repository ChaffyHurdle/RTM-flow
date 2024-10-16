%% Compute kappa_ij
function [kappa_t_elem, bndry_cond] = compute_kappa_ij(t,kappa_elem, ...
                                                grad_lambda,obj_data)

% Various reused objects
dt_vec = diff(obj_data.RTMflow_class.times);
pressure_gradient_tplus1 = obj_data.RTMflow_class.pressure_gradients{t+1};
kappa_t_elem = kappa_elem;
bndry_cond = zeros(1,length(obj_data.mesh_class.nodes));
all_active_elements = obj_data.RTMflow_class.all_active_elements;
centroids = obj_data.mesh_class.centroids;
edge_data = obj_data.RTMflow_class.edge_data;
has_node_i = obj_data.RTMflow_class.Voronoi_mesh_class.has_node_i;
u = obj_data.u;


%% Edges
if t > 10
    moving_boundary_t = obj_data.RTMflow_class.moving_boundary(:,t);
    moving_boundary_inds = find(moving_boundary_t);

    % Find new active elements, assign value based on nearest element at t-1
    new_active_elements = setdiff(edge_data{t}(:,1),edge_data{t+1}(:,1));
    prev_nodes = obj_data.RTMflow_class.moving_boundary(:,t+1);
    prev_nodes_inds = find(prev_nodes);

    prev_elems = zeros(obj_data.mesh_class.num_elements,1);
    for i = 1 : length(prev_nodes_inds)
        candidate = prev_nodes_inds(i);
        for j = 2 : has_node_i(candidate,1)+1
            elem = has_node_i(candidate,j);
            if all_active_elements(elem,t) < 0.5
                prev_elems(elem) = 1;
            end
            %candidate_elem(elem) = 1;
        end
    end

    prev_elems = find(prev_elems);
    nearest_inds = zeros(obj_data.mesh_class.num_elements,1);
    for i = 1:length(new_active_elements)
        elem = new_active_elements(i);
        centroid = centroids(elem,:);
        active_centroids = centroids(prev_elems,:);
        distances = sqrt(sum((active_centroids - centroid).^2, 2));
        [~, idx] = min(distances);
        nearest_ind = prev_elems(idx);
        nearest_inds(elem) = nearest_ind;
        [~,b] = sort(distances);
        kappa_t_elem(elem) = mean( kappa_elem(prev_elems(b(1:3))));
    end    
    
    % Triangles on boundary minus new active elements
    candidate_elem = zeros(obj_data.mesh_class.num_elements,1);
    candidate_elem(edge_data{t}(:,1)) = 1;
    candidate_elem_inds = find(candidate_elem);
    
    for i = 1:length(candidate_elem_inds)
    
        elem_index = candidate_elem_inds(i);
    
        %% extract linearised velocity in element centre
        if ismember(elem_index, new_active_elements)
            nearest_element_ind = nearest_inds(elem);
            grad_p_dot_n = local_flux_tri(edge_data{t}(i,4:5),pressure_gradient_tplus1(nearest_element_ind,:)');
            grad_lambda_dot_n = local_flux_tri(edge_data{t}(i,4:5),grad_lambda(nearest_element_ind,:)');
            dkappa_dt = - exp(u(elem_index)) * grad_p_dot_n * grad_lambda_dot_n;
        else
            grad_p_dot_n = local_flux_tri(edge_data{t}(i,4:5),pressure_gradient_tplus1(elem_index,:)');
            grad_lambda_dot_n = local_flux_tri(edge_data{t}(i,4:5),grad_lambda(elem_index,:)');
            dkappa_dt = - exp(u(elem_index)) * grad_p_dot_n * grad_lambda_dot_n;
        end
        
        kappa_t_elem(elem_index) = kappa_t_elem(elem_index) - dt_vec(t)*dkappa_dt;
    
    end

    for i = 1:length(moving_boundary_inds)
        node_ind = moving_boundary_inds(i);
        rows_involving_node = any(edge_data{t}(:,2:3) == node_ind,2);
        bndry_cond(node_ind) = mean(kappa_t_elem(edge_data{t}(rows_involving_node,1)));
    end
    
end

end


%% Compute v in normal direction
function qn = local_flux_tri(normal_vec,v)

qn = normal_vec*v;

end