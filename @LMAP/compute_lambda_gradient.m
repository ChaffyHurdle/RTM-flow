function grad_lambda_ij = compute_lambda_gradient(obj,lambda_ij)

elements = obj.mesh_class.elements;
num_elements = obj.mesh_class.num_elements;
num_times = length(obj.RTMflow_class.times);
grad_lambda_ij = zeros(num_elements,2,num_times);
shape_gradients = obj.RTMflow_class.pressure_class.shape_fun_gradients;

for i = 1:num_times

    lambda_ij_time = lambda_ij(:,i);

    for j = 1:num_elements 
        
        local_pressure = lambda_ij_time(elements(j,:));
        grad_phi = shape_gradients{j};
    
        grad_lambda_ij(j,:,i) = (grad_phi * local_pressure)';
    end
end

end