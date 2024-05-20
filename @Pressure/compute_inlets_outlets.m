function obj = compute_inlets_outlets(obj)

%% Extract needed mesh information
nodes = obj.mesh_class.nodes;
num_nodes = obj.mesh_class.num_nodes;
boundary_nodes = obj.mesh_class.boundary_nodes;

%% Legacy code inlet/vent flags
obj.inlet_flag = false(num_nodes,1);
obj.vent_flag = false(num_nodes,1);

%% Loop over each boundary node to check for inlet/outlet
for i = 1:length(boundary_nodes)

    node_index = boundary_nodes(i);
    boundary_point = nodes(node_index,:);
    obj.inlet_flag(node_index) = obj.inlet_func(boundary_point);
    obj.vent_flag(node_index) = obj.vent_func(boundary_point); 
end

obj.Dirichlet = obj.inlet_flag;
obj.inlet_pos = nodes(obj.inlet_flag,:);
obj.vent_elem = 0;
obj.vent_idx = find(obj.vent_flag);

%% Find Neumann boundary nodes
obj.Neumann_flag = false(num_nodes,1);
obj.Neumann_flag(boundary_nodes) = true;
obj.Neumann_flag(obj.inlet_flag | obj.vent_flag) = false;
end