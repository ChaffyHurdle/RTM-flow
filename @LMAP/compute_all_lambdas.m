function obj = compute_all_lambdas(obj)

% Various init variables
nobs = obj.RTMflow_class.nobservations;
nsensors = obj.RTMflow_class.nsensors;
i_vec = 1:nsensors;
j_vec = 1:nobs;
[i_vec, j_vec] = meshgrid(i_vec,j_vec);
i_vec = reshape(i_vec,[],1);
j_vec = reshape(j_vec,[],1);
lambda_mat = zeros(length(obj.mesh_class.nodes),length(obj.RTMflow_class.times),nobs*nsensors);
%grad_lambda_mat = cell(1,nobs*nsensors);

% Large objects
stiffness_matrices = obj.RTMflow_class.stiffness_matrices;
active_matrix = obj.RTMflow_class.active_nodes; % nodes in D(t)
dirichlet_matrix = obj.RTMflow_class.Dirichlet_nodes; % inlet, front and beyond front
%shape_cell = obj.RTMflow_class.pressure_class.shape_fun_gradients;

% Compute each adjoint for pressure (along with its gradient)
disp("Entering parallel loop")
parfor k = 1:nobs*nsensors
    i = i_vec(k); j = j_vec(k);

    lambda_ij = obj.compute_lambda_ij(i,j,stiffness_matrices,active_matrix,dirichlet_matrix);
    %grad_lambda_ij = obj.compute_grad_lambda_ij(lambda_ij,shape_cell);
    
    lambda_mat(:,:,k) = lambda_ij;
    %grad_lambda_mat{k} = grad_lambda_ij;

    disp([i,j])
end

% Set data
obj.lambdas = lambda_mat;
%obj.gradient_lambdas = grad_lambda_mat;