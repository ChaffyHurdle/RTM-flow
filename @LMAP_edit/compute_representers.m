function [R,Q] = compute_representers(obj,t)

% Shorthand variables
nodes = obj.mesh_class.nodes;
elements = obj.mesh_class.elements;
times = obj.RTMflow_class.times;
diff_times = diff(times);

num_nodes = obj.mesh_class.num_nodes;
num_elems = obj.mesh_class.num_elements;

stiffness_matrices = obj.RTMflow_class.stiffness_matrices;
active_nodes = obj.RTMflow_class.active_nodes;
Dirichlet_nodes = obj.RTMflow_class.Dirichlet_nodes;
moving_boundary = obj.RTMflow_class.moving_boundary;

grad_pressures = obj.RTMflow_class.pressure_gradients;
active_elements = obj.RTMflow_class.all_active_elements;

phi = obj.physics_class.porosity;
mu = obj.physics_class.viscosity;

% Active sensors
active_sensors = obj.find_active_sensors(t);
starting_triggers = sort(unique(active_sensors(:,6)),'descend');

% Initialise variables
kappas_nodes = zeros(num_nodes,length(active_sensors));
kappas_elems = zeros(num_elems,length(active_sensors));
Q = zeros(num_elems,t*obj.physics_class.nsensors);

% Solve adjoint backwards in time!
k = 1;
active_sensors_t = active_sensors(active_sensors(:,6)>=starting_triggers(k),:);
for t_ind = starting_triggers(k):-1:1
    
    if k < t
        if t_ind <= starting_triggers(k+1)
            k = k + 1;
            active_sensors_t = active_sensors(active_sensors(:,6)>=starting_triggers(k),:);
        end
    end
  
    % Various data from the forward solve
    stiffness_matrix_t = stiffness_matrices{t_ind};
    free = active_nodes(:,t_ind) & ~Dirichlet_nodes(:,t_ind); % nodes in D(t) except front and inlet
    fixed = ~free; % nodes that are not free
    moving_boundary_t = moving_boundary(:,t_ind); % nodes that are active, Dirichlet but not inlet
    
    % Set Dirichlet boundary conditions
    lambda = zeros(length(nodes),size(active_sensors_t,1));
    lambda(moving_boundary_t,:) = 1;
    
    % Compute RHS
    load_vector = obj.compute_load_vec(active_sensors_t,t_ind);

    b_free = load_vector(free,:) - ...
                      stiffness_matrix_t(free,fixed)*lambda(fixed,:);

    % Assemble and solve linear systems
    A_free = stiffness_matrix_t(free,free);
    lambda(free,:) = A_free\b_free;

    % Compute gradients
    grad_lambda = zeros(num_elems,2,size(active_sensors_t,1));
    active_elements_t = active_elements(:,t_ind);
    active_elems_inds = find(active_elements_t);
    shape_fun_grads_t = obj.RTMflow_class.pressure_class.shape_fun_gradients(:,:,active_elems_inds);
    elem_t = obj.mesh_class.elements(active_elems_inds,:);
    for j = 1:size(active_sensors_t,1)
        grad_lambda(active_elems_inds,:,j) = compute_grad_lambda_ij(lambda(:,j), ...
                                                                    shape_fun_grads_t, ...
                                                                    elem_t);
    end

    % subplot(1,2,1)
    % pdeplot(obj.mesh_class.nodes',...
    %     obj.mesh_class.elements', ...
    %     XYData=lambda(:,1), ...
    %     FlowData=grad_lambda(:,:,1), ...
    %     XYStyle='interp',ColorMap="jet",Mesh="off")
    % subplot(1,2,2)
    % pdeplot(obj.mesh_class.nodes',...
    %     obj.mesh_class.elements', ...
    %     XYData=lambda(:,5), ...
    %     FlowData=grad_lambda(:,:,5), ...
    %     XYStyle='interp',ColorMap="jet",Mesh="off")
    % drawnow
    
    grad_pressures_t = grad_pressures{t_ind};
    f_i = active_elements_t .* squeeze(sum(grad_lambda .* grad_pressures_t,2));
    % size(Q(:,active_sensors_t(:,end)))
    % size(f_i)
    Q(:,active_sensors_t(:,end)) = Q(:,active_sensors_t(:,end)) ...
                               + f_i * diff_times(t_ind);
end

Q = - Q .* exp(obj.u');
Q(abs(Q)<1e-10)=0;
R = obj.inverse_class.C0_inv * (Q .* obj.mesh_class.element_areas);

end


function grad_lambda_ij = compute_grad_lambda_ij(lambda_ij,shape_fun_grads,elements)

local_lambdas = lambda_ij(elements);
local_lambdas = permute(local_lambdas, [3, 2, 1]);

part_1 = bsxfun(@times, shape_fun_grads, local_lambdas);
part_2 = sum(part_1,2);
grad_lambda_ij = squeeze(part_2)';

end