function [p,ups] = forward_map(u,params)

options = optimset('Display','off');
upsilon_vec = zeros(1,params.Nt);
pressure_mat = ones(params.Nt,params.Nx)*params.p_0;
p = ones(params.nsensors,params.nobtimes)*params.p_0;
ups = zeros(1,params.nobtimes);
pressure_mat(:,1) = params.p_I;
dt = params.dt;

for i = 2:params.Nt
    opt_func = @(x) x - upsilon_vec(i-1) - (params.p_I-params.p_0)*dt/(F_u(x,u,params)*params.mu*params.phi);
    upsilon = min(fsolve(opt_func,upsilon_vec(i-1),options),params.L);
    pressure = (params.p_I - (params.p_I - params.p_0) * F_u(params.x_locations,u,params)'/F_u(upsilon,u,params)) .* (params.x_locations <= upsilon) ...
        + params.p_0 * (params.x_locations > upsilon);

    upsilon_vec(i) = upsilon;
    pressure_mat(i,:) = pressure;
end

for i = 1:params.nobtimes
    p(:,i) = interp2(params.x_locations, params.t_locations, pressure_mat, params.sensor_locs,ones(1,params.nsensors)*params.ob_times(i), 'linear');
    ups(i) = upsilon_vec(params.obtime_inds(i));
end


end