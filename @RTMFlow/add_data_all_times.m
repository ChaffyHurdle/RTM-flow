function obj = add_data_all_times(obj,it,new_time,new_pressure,...
    new_pressure_gradients,new_flow_rates,new_filling_factors,...
    stiffness_matrix,active_nodes,Dirichlet_nodes,active_elements, new_active_elements,new_filled_volume)

obj.times = [obj.times new_time];
obj.pressures = [obj.pressures new_pressure];
obj.pressure_gradients{it} = new_pressure_gradients;
obj.flow_rates = [obj.flow_rates new_flow_rates];
obj.filling_factors = [obj.filling_factors new_filling_factors];
obj.stiffness_matrices{it} = stiffness_matrix;
obj.active_nodes = [obj.active_nodes active_nodes];
obj.Dirichlet_nodes = [obj.Dirichlet_nodes Dirichlet_nodes];
obj.all_active_elements = [obj.all_active_elements active_elements];
obj.all_new_active_elements = [obj.all_new_active_elements new_active_elements];
obj.new_filled_volumes = [obj.new_filled_volumes new_filled_volume];

end