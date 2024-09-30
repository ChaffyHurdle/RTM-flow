function obj = update_boundary_condition3(obj,t)

% Various reused objects
mesh_class = obj.mesh_class;
dt_vec = diff(obj.RTMflow_u.times);
viscosity = obj.RTMflow_u.physics_class.viscosity;
porosity = obj.RTMflow_u.physics_class.porosity;
pressure_gradient_t = obj.RTMflow_u.pressure_gradients{t};
pressure_gradient_t_minus_1 = obj.RTMflow_u.pressure_gradients{t-1};

%% Edges
if t > 4
    active_nodes = boolean(obj.RTMflow_u.active_nodes(:,t)); % nodes in D(t)
    is_moving_boundary = obj.is_moving_boundary(:,t);
    moving_boundary_inds = find(is_moving_boundary);

    obj = obj.extract_edge_data(t);

    % Find new active elements, assign value based on nearest element at t-1
    new_active_elements = find(obj.RTMflow_u.all_new_active_elements(:,t-1));
    prev_nodes = obj.is_moving_boundary(:,t-1) | obj.RTMflow_u.pressure_class.is_inlet;
    prev_nodes_inds = find(prev_nodes);

    prev_elems = zeros(mesh_class.num_elements,1);
    for i = 1 : length(prev_nodes_inds)
        candidate = prev_nodes_inds(i);
        for j = 2 : obj.RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,1)+1
            elem = obj.RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,j);
            if obj.RTMflow_u.all_active_elements(elem,t-1) >= 0.5
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
        obj.v_h(elem,:) = obj.v_h(nearest_ind,:);
    end

%     if t < 60 && t > 50
% 
%         figure(1)
%         subplot(1,2,1)
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData = obj.RTMflow_u.all_new_active_elements(:,t-1), ...
%             XYStyle='flat',ColorMap="jet",Mesh="on")
%         hold on
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
% 
%         subplot(1,2,2)
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData = nearest_inds, ...
%             XYStyle='flat',ColorMap="jet",Mesh="on")
%         hold on
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
% 
% 
%         figure(2)
%         subplot(1,2,1)
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData=candidate_elem + obj.RTMflow_u.all_new_active_elements(:,t-1), ...
%             XYStyle='flat',ColorMap="jet",Mesh="on");
%         hold on
%         for i = 1:length(edgedata)
%             plot([mesh_class.nodes(edgedata(i,2),1),mesh_class.nodes(edgedata(i,3),1)], ...
%                 [mesh_class.nodes(edgedata(i,2),2),mesh_class.nodes(edgedata(i,3),2)],'w');
%         end
%         hold off
%         hold on
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
% 
%         subplot(1,2,2)
%         all_normals = zeros(mesh_class.num_elements,2);
%         all_normals(edgedata(:,1),:) = normals;
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData=vecnorm(all_normals'), ...
%             XYStyle='flat',ColorMap="jet",Mesh="on");
%         hold on
%         for i = 1:length(edgedata)
%             plot([mesh_class.nodes(edgedata(i,2),1),mesh_class.nodes(edgedata(i,3),1)], ...
%                 [mesh_class.nodes(edgedata(i,2),2),mesh_class.nodes(edgedata(i,3),2)],'w');
%         end
%         for i = 1:length(edgedata)
%             plot([mesh_class.centroids(edgedata(i,1),1),mesh_class.centroids(edgedata(i,1),1)+0.01*normals(i,1)], ...
%                 [mesh_class.centroids(edgedata(i,1),2),mesh_class.centroids(edgedata(i,1),2)+0.01*normals(i,2)],'w');
%         end
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
% 
%         figure(3)
%         subplot(1,2,1)
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData=vecnorm(pressure_gradient_t_minus_1')', ...
%             XYStyle='flat',ColorMap="jet",Mesh="on");
%         hold on
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
%         subplot(1,2,2)
%         pdeplot(mesh_class.nodes',mesh_class.elements',XYData=vecnorm(obj.grad_p_tilde')', ...
%             XYStyle='flat',ColorMap="jet",Mesh="on");
%         hold on
%         plot(mesh_class.nodes(moving_boundary_inds,1),mesh_class.nodes(moving_boundary_inds,2),'wo')
%         hold off
%         input('')
%     end
    
    
    % Triangles on boundary minus new active elements
    candidate_elem = zeros(mesh_class.num_elements,1);
    candidate_elem(obj.edge_data(:,1)) = 1;
    candidate_elem_inds = find(candidate_elem);

    bnd_vals = zeros(mesh_class.num_elements,1);
    
    for i = 1:length(candidate_elem_inds)
    
        elem_index = candidate_elem_inds(i);
    
        %% extract linearised velocity in element centre
        if ismember(elem_index, new_active_elements)
            nearest_element_ind = nearest_inds(elem);
            dvh_dt = - (1/(viscosity*porosity)) * ( exp(obj.u(elem_index)) * obj.h(elem_index) * pressure_gradient_t_minus_1(nearest_element_ind,:)' ...
                              + exp(obj.u(elem_index)) * obj.grad_p_tilde(nearest_element_ind,:)');
        else
            dvh_dt = - (1/(viscosity*porosity)) * ( exp(obj.u(elem_index)) * obj.h(elem_index) * pressure_gradient_t_minus_1(elem_index,:)' ...
                              + exp(obj.u(elem_index)) * obj.grad_p_tilde(elem_index,:)' );
        end
        
        obj.v_h(elem_index,:) = obj.v_h(elem_index,:) + dt_vec(t-1)*dvh_dt';
    
        %% Compute transfer into each CV
        local_linearised_flow_rate = local_flux_tri(obj.edge_data(i,4:5),obj.v_h(elem_index,:)');
        local_pressure_grad = local_flux_tri(obj.edge_data(i,4:5),pressure_gradient_t(elem_index,:)');
    
        %% Need to track v_h as it depends on time
        bnd_vals(elem_index) = - local_linearised_flow_rate*local_pressure_grad;
    
    end
    %obj.v_h(obj.RTMflow_u.pressure_class.is_inlet,:) = 0;

    for i = 1:length(moving_boundary_inds)
        node_ind = moving_boundary_inds(i);
        rows_involving_node = any(obj.edge_data(:,2:3) == node_ind,2);
        obj.bndry_cond(node_ind) = mean(bnd_vals(obj.edge_data(rows_involving_node,1)));
    end
    
%     times_u_plus_h = obj.RTMflow_u_plus_h.times;
%     nearest_t_u_plus_h = find(times_u_plus_h > obj.RTMflow_u_plus_h.times(t), 1);
%     obj.bndry_cond(obj.is_moving_boundary(:,t) & obj.is_moving_boundary_u_plus_h(:,nearest_t_u_plus_h-1)) = 0;
    %obj.bndry_cond(obj.is_moving_boundary(:,t) & ~obj.RTMflow_u_plus_h.active_nodes(:,nearest_t_u_plus_h-1)) = 0;
    
%     figure(1)
%     pdeplot(mesh_class.nodes',...
%            mesh_class.elements', ...
%            XYData=bnd_vals,XYStyle='flat',ColorMap="jet",Mesh="on")
    
end

end

function qn = local_flux_tri(normal_vec,v)

qn = normal_vec*v;

end