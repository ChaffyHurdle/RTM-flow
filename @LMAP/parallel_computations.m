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

addAttachedFiles(gcp,["compute_kappa_ij.m","compute_representer_ij.m"])
parfor k = 1 : nobs * nsensors
    i = obj_data.i_vec(k);
    j = obj_data.j_vec(k);
    %disp([i,j])

    t_j = obj_data.physics_class.observation_times(j);
    x_i_sensor_elem = obj_data.RTMflow_class.sensor_element_inds(i);
    time_index = find(obj_data.RTMflow_class.times > t_j,1);
    sensor_active = ismember(x_i_sensor_elem,find(obj_data.RTMflow_class.all_active_elements(:,time_index)));

    if sensor_active
       % do stuff
       [lambda_ij,grad_lambda_ij] = compute_lambda_ij(i,j,obj_data);
       lambda_mat(:,:,k) = lambda_ij;
        
       Q_ij = compute_representer_ij(grad_lambda_ij,obj_data);
       R_ij = obj_data.inverse_class.C0_inv * (Q_ij' .* obj_data.mesh_class.element_areas);
       Q_mat(:,k) = Q_ij;
       R_mat(:,k) = R_ij;
       
       %ptilde_ij = compute_ptilde_ij(R_ij,obj_data);
       %p_tilde_mat(:,:,k) = ptilde_ij;
    end
end
%R_0 = obj.u0 - obj.u;
%obj.p_tilde0 = compute_ptilde_ij(R_0,obj_data);
obj.lambdas = lambda_mat;
obj.Q = Q_mat;
obj.R = R_mat;
%obj.p_tildes = p_tilde_mat;
end