function [lambda_ij,grad_lambda_ij] = compute_lambda_ij(i,j,obj_data)

% Locate (x_i,t_j)
t_j = obj_data.physics_class.observation_times(j);
x_i = obj_data.physics_class.sensor_locs(i,:);
x_i_sensor_elem = obj_data.RTMflow_class.sensor_element_inds(i);

% Shorthand variables
nodes = obj_data.mesh_class.nodes;
num_nodes = obj_data.mesh_class.num_nodes;
elements = obj_data.mesh_class.elements;
num_elems = obj_data.mesh_class.num_elements;
times = obj_data.RTMflow_class.times;
stiffness_matrices = obj_data.RTMflow_class.stiffness_matrices;
active_nodes = obj_data.RTMflow_class.active_nodes;
Dirichlet_nodes = obj_data.RTMflow_class.Dirichlet_nodes;
moving_boundary = obj_data.RTMflow_class.moving_boundary;
phi = obj_data.physics_class.porosity;
mu = obj_data.physics_class.viscosity;

% Initialise variables
lambda_ij = zeros(num_nodes,length(times));
grad_lambda_ij = zeros(num_elems,2,length(times));
kappa_elem_tplus1 = zeros(1,length(elements));
grad_lambda = zeros(length(elements),2);
[~,closest_time] = sort(abs(obj_data.RTMflow_class.times-t_j));
closest_time = closest_time(1);
del_t = max(obj_data.RTMflow_class.times(closest_time+2) - obj_data.RTMflow_class.times(closest_time),...
            obj_data.RTMflow_class.times(closest_time) - obj_data.RTMflow_class.times(closest_time-2));
del_t = del_t^2;
start_index = min(closest_time+100,length(times)-1);
%del_t = 0.002/200;

% Compute aspect of load vector from \delta(x_i-x) contribution.
nodes_surrounding_xi = nodes(elements(x_i_sensor_elem,:),:);
[lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding_xi, x_i);
load_vector_x = zeros(num_nodes,1);
load_vector_x(elements(x_i_sensor_elem,1)) = lambda1;
load_vector_x(elements(x_i_sensor_elem,2)) = lambda2;
load_vector_x(elements(x_i_sensor_elem,3)) = lambda3;

% Solve adjoint backwards in time!
for t = start_index:-1:1
  
    % Various data from the forward solve
    stiffness_matrix_t = stiffness_matrices{t};
    free = active_nodes(:,t) & ~Dirichlet_nodes(:,t); % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free
    moving_boundary_t = moving_boundary(:,t); % nodes that are active, Dirichlet but not inlet
    
    % Set Dirichlet boundary conditions
    lambda = zeros(length(nodes),1);
    [kappa_elem, bndry_cond] = compute_kappa_ij(t,kappa_elem_tplus1,grad_lambda,obj_data);
    lambda(moving_boundary_t) = bndry_cond(moving_boundary_t)/(mu*phi);
    
    load_vector_t = obj_data.delta_t(t_j,times(t),del_t);
    load_vector = load_vector_t*load_vector_x;

    b_free = load_vector(free) - ...
                      stiffness_matrix_t(free,fixed)*lambda(fixed);
    A_free = stiffness_matrix_t(free,free);

    % solve matrix equation for lambda_ij
    lambda(free) = A_free\b_free;
    lambda_ij(:,t) = lambda;
    kappa_elem_tplus1 = kappa_elem;

    % Compute gradient lambda_ij
    grad_lambda = compute_grad_lambda_ij(lambda,obj_data);
    grad_lambda_ij(:,:,t) = grad_lambda;
end
end


%% Compute barycentric coords
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p(1) - v3(1)) + (v3(1) - v2(1))*(p(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p(1) - v3(1)) + (v1(1) - v3(1))*(p(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end

function grad_lambda_ij = compute_grad_lambda_ij(lambda_ij,obj_data)

elements = obj_data.mesh_class.elements;
shape_fun_grads = obj_data.RTMflow_class.pressure_class.shape_fun_gradients;

local_lambdas = lambda_ij(elements);
grad_lambda_ij = squeeze(sum(bsxfun(@times, shape_fun_grads, permute(local_lambdas, [3, 2, 1])), 2))';

end