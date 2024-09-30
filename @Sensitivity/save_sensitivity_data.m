function obj = save_sensitivity_data(obj)

obj.p_tildes = obj.p_tilde(:,obj.time_inds_u);
obj.is_moving_boundary_ob_times = obj.is_moving_boundary(:,obj.time_inds_u);
obj.is_moving_boundary_ob_times_u_plus_h = obj.is_moving_boundary_u_plus_h(:,obj.time_inds_u_plus_h);

end