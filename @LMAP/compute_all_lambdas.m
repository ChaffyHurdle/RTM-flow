function obj = compute_all_lambdas(obj)

nobs = obj.RTMflow_class.nobservations;
nsensors = obj.RTMflow_class.nsensors;
i_vec = 1:nsensors;
j_vec = 1:nobs;
[i_vec, j_vec] = meshgrid(i_vec,j_vec);
i_vec = reshape(i_vec,[],1);
j_vec = reshape(j_vec,[],1);
lambda_cell = zeros(obj.mesh_class.num_nodes,length(obj.RTMflow_class.times),nobs*nsensors);

parfor k = 1:nobs*nsensors
    
    i = i_vec(k);
    j = j_vec(k);
    lambda_ij = obj.compute_lambda_ij(i,j);
    disp([i,j,size(lambda_ij)])
    lambda_cell(:,:,k) = lambda_ij;
end

obj.lambdas = lambda_cell;