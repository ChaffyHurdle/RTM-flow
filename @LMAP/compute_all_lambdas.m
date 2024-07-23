function obj = compute_all_lambdas(obj)

nobs = obj.RTMflow_class.nobservations;
nsensors = obj.RTMflow_class.nsensors;
i_vec = 1:nsensors;
j_vec = 1:nobs;
[i_vec, j_vec] = meshgrid(i_vec,j_vec);
i_vec = reshape(i_vec,[],1);
j_vec = reshape(j_vec,[],1);
lambda_mat = cell(1,nobs*nsensors);
grad_lambda_mat = cell(1,nobs*nsensors);

parfor k = 1:nobs*nsensors
    
    i = i_vec(k);
    j = j_vec(k);
    lambda_ij = obj.compute_lambda_ij(i,j);
    grad_lambda_ij = obj.compute_lambda_gradient(lambda_ij);
    disp([i,j,size(lambda_ij)])
    lambda_mat{k} = lambda_ij;
    grad_lambda_mat{k} = grad_lambda_ij;
end

obj.lambdas = lambda_mat;
obj.gradient_lambdas = grad_lambda_mat;