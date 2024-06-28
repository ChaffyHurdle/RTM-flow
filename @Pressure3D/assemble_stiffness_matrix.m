function obj = assemble_stiffness_matrix(obj)
%% This function constructs the stiffness matrix required to solve for the 
%% pressure field.

new_elements = find(obj.new_active_elements);

K = obj.physics_class.permeability;

for i = 1 : length(new_elements)

    %% element geometry
    nodes = obj.mesh_class.elements(new_elements(i),:);
    centroid = obj.mesh_class.centroids(new_elements(i),:);
    measure_element = obj.mesh_class.element_measures(new_elements(i));
    
    %% local FEM stiffness matrix
    grad_phi = obj.shape_fun_gradients{new_elements(i)};
    A_local = (K(centroid)*grad_phi)'*grad_phi*measure_element;
    
    %% local to global mapping
    obj.stiffness_matrix(nodes,nodes) = ...
                               obj.stiffness_matrix(nodes,nodes) + A_local;
end
end
