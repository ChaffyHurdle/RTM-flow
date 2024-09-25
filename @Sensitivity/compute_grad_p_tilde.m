function obj = compute_grad_p_tilde(obj,p_tilde)

elements = obj.mesh_class.elements;
num_elements = obj.mesh_class.num_elements;

grad_p_tilde = zeros(num_elements,2);
for i = 1:num_elements 
    local_p_tilde = p_tilde(elements(i,:));
    grad_phi = obj.RTMflow_u.pressure_class.shape_fun_gradients{i};
    grad_p_tilde(i,:) = (grad_phi * local_p_tilde)';
end
obj.grad_p_tilde = grad_p_tilde;
end