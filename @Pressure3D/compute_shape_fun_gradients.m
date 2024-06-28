function obj = compute_shape_fun_gradients(obj)

num_elements = obj.mesh_class.num_elements;

obj.shape_fun_gradients = cell(num_elements,1);

for i = 1:num_elements

    element_nodes = obj.mesh_class.elements(i,:);

    x = obj.mesh_class.nodes(element_nodes,1); 
    y = obj.mesh_class.nodes(element_nodes,2);
    z = obj.mesh_class.nodes(element_nodes,3);
    
    J = [x(2)-x(1) x(3)-x(1) x(4)-x(1);
         y(2)-y(1) y(3)-y(1) y(4)-y(1);
         z(2)-z(1) z(3)-z(1) z(4)-z(1)];

    grad_ref = [-1 1 0 0;
                -1 0 1 0;
                -1 0 0 1];
    
    grad_phi = (J')\grad_ref;

    obj.shape_fun_gradients{i} = grad_phi;
end

end