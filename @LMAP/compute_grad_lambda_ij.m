function grad_lambda_ij = compute_grad_lambda_ij(obj,lambda_ij,shape_cell)

elements = obj.mesh_class.elements;
num_elements = obj.mesh_class.num_elements;
num_times = length(obj.RTMflow_class.times);
grad_lambda_ij = zeros(num_elements,2,num_times);

for i = 1:num_times

    lambda_ij_time = lambda_ij(:,i);

    for j = 1:num_elements 
        
        local_lambda = lambda_ij_time(elements(j,:));
        grad_phi = shape_cell.Values{j};
    
        grad_lambda_ij(j,:,i) = (grad_phi * local_lambda)';
    end
end

end