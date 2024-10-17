function obj = parallel_computations(obj)

% Shorthand variables
num_elems = obj.mesh_class.num_elements;
num_nodes = obj.mesh_class.num_nodes;
num_times = length(obj.RTMflow_class.times);
nobs = length(obj.physics_class.observation_times);
nsensors = obj.physics_class.nsensors;

% Initialise outputs
lambda_mat = zeros(num_nodes,num_times, nobs * nsensors);
R_mat = zeros(num_elems, nobs * nsensors);
Q_mat = zeros(num_elems, nobs * nsensors);

% Send data to all workers
obj_data = obj;

% Find all non-trivial indices
active_sensors = zeros(1,nobs*nsensors);
for k = 1:nobs*nsensors
    i = obj_data.i_vec(k);
    j = obj_data.j_vec(k);
    t_j = obj_data.physics_class.observation_times(j);
    x_i_sensor_elem = obj_data.RTMflow_class.sensor_element_inds(i);
    time_index = find(obj_data.RTMflow_class.times > t_j,1);
    active_sensors(k) = ismember(x_i_sensor_elem,find(obj_data.RTMflow_class.all_active_elements(:,time_index)));
end

active_inds = find(active_sensors);
active_times = obj_data.j_vec(active_inds);
[~,sorted_active_times_inds] = sort(active_times,'descend');
sorted_inds = active_inds(sorted_active_times_inds);
num_taken = cpu_scheduler(obj_data,active_sensors);
cum_sum_num_taken = cumsum(num_taken);

% Creat sub-matrices of active indices
active_lambda_mat = zeros(num_nodes,num_times, length(active_inds));
active_R_mat = zeros(num_elems, length(active_inds));
active_Q_mat = zeros(num_elems, length(active_inds));

for i = 1:length(num_taken)-1
    % disp("Starting parfor " + num2str(i))
    parfor k = cum_sum_num_taken(i)+1 : cum_sum_num_taken(i+1)
        idx = sorted_inds(k);
        i = obj_data.i_vec(idx);
        j = obj_data.j_vec(idx);
        [lambda_ij,grad_lambda_ij] = compute_lambda_ij(i,j,obj_data);
        active_lambda_mat(:,:,k) = lambda_ij;
        
        Q_ij = compute_representer_ij(grad_lambda_ij,obj_data);
        R_ij = obj_data.inverse_class.C0_inv * (Q_ij' .* obj_data.mesh_class.element_areas);
        active_Q_mat(:,k) = Q_ij;
        active_R_mat(:,k) = R_ij;
    end
end
% disp("Finished parfor")

% Add sub-matrix parts to overall matrix
for k = 1 : length(sorted_inds)
    lambda_mat(:,:,sorted_inds(k)) = active_lambda_mat(:,:,k);
    Q_mat(:,sorted_inds(k)) = active_Q_mat(:,k);
    R_mat(:,sorted_inds(k)) = active_R_mat(:,k);
end

obj.lambdas = lambda_mat;
obj.Q = Q_mat;
obj.R = R_mat;
%obj.p_tildes = p_tilde_mat;
end