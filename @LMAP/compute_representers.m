function obj = compute_representers(obj)

% Various init variables
nobs = obj.physics_class.nobservations;
nsensors = obj.physics_class.nsensors;
representers = zeros(obj.mesh_class.elements,nsensors*nobs);

parfor k = 1:nobs * nsensors
    
    R_i = compute_representer_i(k);
    representers(:,k) = R_i;
end
obj.representers = representers;

end


function R_i = compute_representer_i(i)
    R_i = i;
end