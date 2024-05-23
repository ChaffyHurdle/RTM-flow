function obj = compute_volume_measures(obj,thickness)
            
obj.volume_measures = zeros(obj.mesh_class.num_nodes,1);

for i = 1:obj.mesh_class.num_elements
    
    element_area = obj.mesh_class.element_areas(i);
    control_volumes = obj.mesh_class.elements(i,:);
    sub_volume_measures = element_area/3 * thickness *ones(3,1);
    obj.volume_measures(control_volumes) = obj.volume_measures(control_volumes) + sub_volume_measures;
end
end