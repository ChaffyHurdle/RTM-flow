function [u_maps,std_maps,iteration_vec,timer_vec] = run_LMAP_sequential(u0,C,params,experiment,plotting)

C_minus_half = inv(sqrtm(C));
tol_J = 0.01;
tol_U = 0.01;
u_maps = zeros(params.Nx,params.nobtimes);
std_maps = zeros(params.Nx,params.nobtimes);
iteration_vec = zeros(1,params.nobtimes);
timer_vec = zeros(1,params.nobtimes);
hex = "#d3d3d3";

for t = 1:params.nobtimes
    alpha = 1e4;
    breaker1 = 0;
    breaker2 = 0;
    breaker3 = 0;
    iterations = 0;
    patience = 0;

    Sigma_t = reshape(experiment.Sigma(:,1:t),[],1);
    d_t = reshape(experiment.d(:,1:t),[],1);
    Sigma_minus_half = diag(1./sqrt(Sigma_t));
    Sigma_t = diag(Sigma_t);

    u = u0;
    [p,ups] = forward_map(u,params);
    data_misfit = norm(Sigma_minus_half * (d_t - reshape(p(:,1:t),[],1)))^2;
    distance_from_u0 = norm(C_minus_half * (u - u0)')^2;
    J = data_misfit + distance_from_u0;

    tic;
    for k = 1:100
        [boldR,curlyR,c] = compute_representers(u,ups,C,params,t);
    
        h = (u0-u)/(1+alpha) + (boldR * ((curlyR + (1+alpha)*Sigma_t)\(d_t - reshape(p(:,1:t),[],1) - c/(1+alpha))))';
    
        u_cand = u + h;
    
        [p_cand,ups_cand] = forward_map(u_cand,params);
        data_misfit = norm(Sigma_minus_half * (d_t - reshape(p_cand(:,1:t),[],1)))^2;
        distance_from_u0 = norm(C_minus_half * (u_cand - u0)')^2;
        J_cand = data_misfit + distance_from_u0;
    
        if J_cand < J
            breaker1 = (max(abs(u-u_cand)./abs(u)) < tol_U);
            breaker2 = (J - J_cand < tol_J);
            u = u_cand;
            J = J_cand;
            ups = ups_cand;
            p = p_cand;
            iterations = iterations + 1;
            alpha = alpha/2;
            patience = 0;
        else
            alpha = alpha*2;
            patience = patience + 1;
            breaker3 = (patience == 5);
        end
    
        if breaker1 || breaker2 || breaker3
            u_map = u;
            C_map = C - boldR * inv(curlyR + Sigma_t) * boldR';
            timer = toc

            u_maps(:,t) = u_map;
            std_maps(:,t) = sqrt(diag(C_map));
            iteration_vec(t) = iterations;
            timer_vec(t) = timer;
            break
        end
    end

end

if plotting
    
    figure(6)
    subplot(2,3,1)
    upper50 = u0 + 0.674*sqrt(diag(C))';
    lower50 = u0 - 0.674*sqrt(diag(C))';
    upper95 = u0 + 1.96*sqrt(diag(C))';
    lower95 = u0 - 1.96*sqrt(diag(C))';
    fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none')
    hold on
    plot(params.x_locations,u0,"k--")
    plot(params.x_locations,upper95,"k")
    plot(params.x_locations,lower95,"k")
    plot(experiment.params_fwd.x_locations,experiment.u_true,"r")
    hold off
    xlim([0,1])
    ylim([-2,2])
    xlabel('$$x$$','interpreter','latex')
    ylabel('$$u(x)$$','interpreter','latex')
    title("$$\mu_0$$",'Interpreter','latex')

    for i = 1:params.nobtimes
        subplot(2,3,i+1)
        upper50 = u_maps(:,i)' + 0.674*std_maps(:,i)';
        lower50 = u_maps(:,i)' - 0.674*std_maps(:,i)';
        upper95 = u_maps(:,i)' + 1.96*std_maps(:,i)';
        lower95 = u_maps(:,i)' - 1.96*std_maps(:,i)';
        fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none')
        hold on
        plot(params.x_locations,u_maps(:,i),"k--")
        plot(params.x_locations,upper95,"k")
        plot(params.x_locations,lower95,"k")
        plot(experiment.params_fwd.x_locations,experiment.u_true,"r")
        xline(experiment.ups_true(i),"r--")
        hold off
        xlim([0,1])
        ylim([-2,2])
        xlabel('$$x$$','interpreter','latex')
        ylabel('$$u(x)$$','interpreter','latex')
        title(sprintf('$\\mu_{%d}$',i),'Interpreter','latex')
    end


end



end