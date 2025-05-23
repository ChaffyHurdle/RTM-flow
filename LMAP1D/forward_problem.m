function [pressure_mat,upsilon_vec] = forward_problem(u,params,plotting)

options = optimset('Display','off');
upsilon_vec = zeros(1,params.Nt);
pressure_mat = ones(params.Nt,params.Nx)*params.p_0;
pressure_mat(:,1) = params.p_I;
dt = params.dt;

for i = 2:params.Nt
    opt_func = @(x) x - upsilon_vec(i-1) - dt/(F_u(x,u,params)*params.mu*params.phi);
    upsilon = min(fsolve(opt_func,upsilon_vec(i-1),options),params.L);
    pressure = (params.p_I - (params.p_I - params.p_0) * F_u(params.x_locations,u,params)'/F_u(upsilon,u,params)) .* (params.x_locations <= upsilon) ...
        + params.p_0 * (params.x_locations > upsilon);

    upsilon_vec(i) = upsilon;
    pressure_mat(i,:) = pressure;
end

if plotting
    figure
    subplot(1,2,1)
    plot(params.x_locations,u);
    ylim([-1.5,1.5]);
    xlim([0,params.L]);
    xlabel("$x$",'interpreter','latex')
    ylabel("$u(x)$",'interpreter','latex')
    title("$u(x)$",'interpreter','latex')

    subplot(1,2,2)
    pcolor(params.x_locations,params.t_locations,pressure_mat);
    xlim([0,params.L]);
    ylim([0,params.T]);
    shading interp;
    grid off;
    colormap jet;
    colorbar
    xlabel("$x$",'interpreter','latex')
    ylabel("$t$",'interpreter','latex')
    title("$p(x,t)$",'interpreter','latex')
    [sensor_grid_x,sensor_grid_y] = meshgrid(params.sensor_locs,params.ob_times);
    hold on
    scatter(sensor_grid_x,sensor_grid_y,'k','filled')
    hold off
end

end