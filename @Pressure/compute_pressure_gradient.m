function obj = compute_pressure_gradient(obj)

elements = obj.mesh_class.elements;
%num_elements = obj.mesh_class.num_elements;

local_pressures = obj.pressure(elements);
obj.pressure_gradient = squeeze(sum(bsxfun(@times, obj.shape_fun_gradients, permute(local_pressures, [3, 2, 1])), 2))';

% for i = 1:num_elements 
% 
%     local_pressure = obj.pressure(elements(i,:));
%     grad_phi = obj.shape_fun_gradients(:,:,i);
% 
%     obj.pressure_gradient(i,:) = (grad_phi * local_pressure)';
% end

end