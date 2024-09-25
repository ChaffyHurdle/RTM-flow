function obj = update_boundary_condition2(obj,t)

% Various variable sets for brevity
mesh_class = obj.mesh_class;
dt_vec = diff(obj.RTMflow_u.times);
viscosity = obj.RTMflow_u.physics_class.viscosity;
porosity = obj.RTMflow_u.physics_class.porosity;
volume_measures = obj.RTMflow_u.Voronoi_mesh_class.volume_measures;
normal_vec = obj.RTMflow_u.Voronoi_mesh_class.volume_outflow_vectors;
pressure_gradient_t = obj.RTMflow_u.pressure_gradients{t};
pressure_gradient_t_minus_1 = obj.RTMflow_u.pressure_gradients{t-1};

% Find nodes associated with inlet, and the boundary at t
is_moving_boundary = obj.is_moving_boundary(:,t);
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
        if obj.RTMflow_u.all_active_elements(elem,t) >= 0.5
            candidate_elem(elem) = 1;
        end
        %candidate_elem(elem) = 1;
    end
end

new_active_elements = find(obj.RTMflow_u.all_new_active_elements(:,t));
repeat_elems = candidate_elem & (~obj.RTMflow_u.all_new_active_elements(:,t));
repeat_elems_ind = find(repeat_elems);

% figure(1)
% subplot(1,3,1)
% pdeplot(mesh_class.nodes',mesh_class.elements',XYData=candidate_elem,XYStyle='flat',ColorMap="jet",Mesh="on")
% subplot(1,3,2)
% pdeplot(mesh_class.nodes',mesh_class.elements',XYData=obj.RTMflow_u.all_new_active_elements(:,t),XYStyle='flat',ColorMap="jet",Mesh="on")
% subplot(1,3,3)
% pdeplot(mesh_class.nodes',mesh_class.elements',XYData=double(candidate_elem & (~obj.RTMflow_u.all_new_active_elements(:,t))),XYStyle='flat',ColorMap="jet",Mesh="on")
% if t > 50 && t < 75
%     input('')
% end

for i = 1:length(new_active_elements)
    elem = new_active_elements(i);
    centroid = mesh_class.centroids(elem,:);
    active_centroids = mesh_class.centroids(repeat_elems,:);
    distances = sqrt(sum((active_centroids - centroid).^2, 2));
    [~, idx] = min(distances);
    nearest_ind = repeat_elems_ind(idx);
    obj.v_h(elem,:) = obj.v_h(nearest_ind,:);
end

candidate_elem = find(candidate_elem);

% Set pressure gradient to 0
grap_p_dot_n = zeros(mesh_class.num_nodes,1);
vh_dot_n = zeros(mesh_class.num_nodes,1);

for i = 1:length(candidate_elem)

    elem_index = candidate_elem(i);

    %% extract local element properties
    element = mesh_class.elements(elem_index,:);

    %% extract linearised velocity in element centre
    dvh_dt = - (1/(viscosity*porosity)) * ( exp(obj.u(elem_index)) * obj.h(elem_index) * pressure_gradient_t_minus_1(elem_index,:)' ...
                      + exp(obj.u(elem_index)) * obj.grad_p_tilde(elem_index,:)' );
    obj.v_h(elem_index,:) = obj.v_h(elem_index,:) + dt_vec(t-1)*dvh_dt';

    %% Convert normal
    normals = normal_vec(elem_index,:);
    normals(1:2) = normals(1:2) / norm(normals(1:2));  
    normals(3:4) = normals(3:4) / norm(normals(3:4));
    normals(5:6) = normals(5:6) / norm(normals(5:6));

    %% Compute transfer into each CV
    local_linearised_flow_rate = local_flux_tri(normal_vec(elem_index,:),obj.v_h(elem_index,:)');
    local_pressure_grad = local_flux_tri(normals,pressure_gradient_t(elem_index,:)');

    %% Need to track v_h as it depends on time
    vh_dot_n(element) = vh_dot_n(element) + local_linearised_flow_rate;
    grap_p_dot_n(element) = grap_p_dot_n(element) + local_pressure_grad;

end
%obj.v_h(inlet_nodes,:) = 0;
obj.bndry_cond = - vh_dot_n .* grap_p_dot_n;
% figure(1)
% subplot(1,4,1)
% pdeplot(mesh_class.nodes',mesh_class.elements',...
%     XYData=(obj.RTMflow_u_plus_h.pressures(:,t)-obj.RTMflow_u.pressures(:,t)),XYStyle='flat',ColorMap="jet",Mesh="on")
% caxis([-0.1,0.1])
% subplot(1,4,2)
% pdeplot(mesh_class.nodes',mesh_class.elements',...
%     XYData=vh_dot_n.*is_moving_boundary,XYStyle='flat',ColorMap="jet",Mesh="on")
% subplot(1,4,3)
% pdeplot(mesh_class.nodes',mesh_class.elements',...
%     XYData=-grap_p_dot_n.*is_moving_boundary,XYStyle='flat',ColorMap="jet",Mesh="on")
% subplot(1,4,4)
% pdeplot(mesh_class.nodes',mesh_class.elements',...
%     XYData=-vh_dot_n.*grap_p_dot_n.*is_moving_boundary,XYStyle='flat',ColorMap="jet",Mesh="on")


end

%% function to compute outflow of standard elements
function qn = local_flux_tri(normal_vec,v)
qn = zeros(3,1);

qn(1) = normal_vec(:,[1 2])*v;
qn(2) = normal_vec(:,[3 4])*v;
qn(3) = normal_vec(:,[5 6])*v;

end