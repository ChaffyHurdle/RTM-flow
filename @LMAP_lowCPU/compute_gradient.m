function grad_lambda_ij = compute_gradient(obj,f)

elements = obj.mesh_class.elements;
shape_fun_grads = obj.RTMflow_class.pressure_class.shape_fun_gradients;

local_f = f(elements);
grad_lambda_ij = squeeze(sum(bsxfun(@times, shape_fun_grads, permute(local_f, [3, 2, 1])), 2))';

end