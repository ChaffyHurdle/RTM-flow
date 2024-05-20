function obj = assemble_stiffness_matrix(obj, newElement)
%% This function constructs the stiffness matrix required to solve for the 
%% pressure field.

new_elements = find(newElement);

K = obj.darcxy_class.permeability;

for i = 1 : length(new_elements)

    %% element geometry
    tri_nodes = obj.mesh_class.elements(new_elements(i),:);
    centroid = obj.mesh_class.centroids(i,:);
    area = obj.mesh_class.element_areas(i);
    
    %% local FEM stiffness matrix
    grad_phi = obj.shape_fun_gradients{i};
    A_local = (K(centroid)*grad_phi)'*grad_phi*area;
    
    %% local to global mapping
    obj.stiffness_matrix(tri_nodes,tri_nodes) = ...
                       obj.stiffness_matrix(tri_nodes,tri_nodes) + A_local;
end
end
