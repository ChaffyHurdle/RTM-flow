function obj = find_sensor_locs_mesh(obj)

% Find nearest sensor locations on computational mesh (indices)
obj.sensor_inds_on_mesh = dsearchn(obj.Delaunay_mesh_class.nodes,obj.physics_class.sensor_locs);

% Locate sensors on computational mesh
obj.sensor_locs_on_mesh = obj.Delaunay_mesh_class.nodes(obj.sensor_inds_on_mesh,:);

end