function obj = add_data_all_times(obj,it,new_time,new_pressure,...
    new_pressure_gradients,new_flow_rates,new_filling_factors,...
    stiffness_matrix,active_nodes,Dirichlet_nodes,active_elements,...
    new_active_elements,new_filled_volume,edge_data)

obj.times(it) = new_time;
obj.pressures(:,it) = new_pressure;
obj.pressure_gradients{it} = new_pressure_gradients;
obj.flow_rates(:,it) = new_flow_rates;
obj.filling_factors(:,it) = new_filling_factors;
obj.stiffness_matrices{it} = stiffness_matrix;
obj.active_nodes(:,it) = active_nodes;
obj.Dirichlet_nodes(:,it) = Dirichlet_nodes;
obj.all_active_elements(:,it) = active_elements;
obj.new_filled_volumes = [obj.new_filled_volumes new_filled_volume];
obj.edge_data{it} = edge_data;

if ~isempty(new_active_elements)
    obj.all_new_active_elements(:,it) = new_active_elements;
end
end