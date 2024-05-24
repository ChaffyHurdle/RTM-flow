function obj = compute_velocity(obj,pressure_class)
            
%% pressure gradients
grad_p = pressure_class.pressure_gradient;

centroids = pressure_class.mesh_class.centroids;

%% Darcy's law applied
K = obj.darcy_class.permeability;
phi = obj.darcy_class.porosity;
mu = obj.darcy_class.viscosity;

for i=1:size(centroids,1)

    obj.velocity(i,:) = -K(centroids(i,:))*grad_p(i,:)'/(phi*mu);
end

end