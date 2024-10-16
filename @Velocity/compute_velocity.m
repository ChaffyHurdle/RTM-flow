function obj = compute_velocity(obj,pressure_class)
            
%% pressure gradients
grad_p = pressure_class.pressure_gradient;

%% Darcy's law applied
K = obj.physics_class.permeability;
phi = obj.physics_class.porosity;
mu = obj.physics_class.viscosity;
obj.velocity = -K.*grad_p/(phi*mu);

% for i=1:size(centroids,1)
%     obj.velocity(i,:) = -K(i)*grad_p(i,:)'/(phi*mu);
% end

end