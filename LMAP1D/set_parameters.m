function params = set_parameters(L,p_I,p_0,mu,phi,T,...
                                 var_matern,l,nu,...
                                 Nx,Nt,...,
                                 nsensors,ob_times)

params.L = L;
params.p_I = p_I;
params.p_0 = p_0;
params.mu = mu;
params.phi = phi;
params.T = T;

params.var_matern = var_matern;
params.l = l;
params.nu = nu;

params.Nx = Nx;
params.Nt = Nt;

params.t_locations = linspace(0,T,Nt);
params.x_locations = linspace(0,L,Nx);
diff_t = diff(params.t_locations);
diff_x = diff(params.x_locations);
params.dt = diff_t(1);
params.dx = diff_x(1);

params.sensor_locs = linspace(0.05,L-0.05,nsensors);
params.nsensors = nsensors;
params.ob_times = ob_times;
params.nobtimes = length(ob_times);

params.obtime_inds = zeros(1,length(ob_times));
for i = 1:length(ob_times)
    [~, idx] = min(abs(params.t_locations - ob_times(i)));
    params.obtime_inds(i) = idx;
end

params.sensor_inds = zeros(1,nsensors);
for i = 1:nsensors
    [~, idx] = min(abs(params.x_locations - params.sensor_locs(i)));
    params.sensor_inds(i) = idx;
end

end