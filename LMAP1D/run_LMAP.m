function [u_map,C_map,timer,iterations] = run_LMAP(u0,C,params,experiment,plotting)

C_minus_half = inv(sqrtm(C));
alpha = 1e4;
tol_J = 0.1;
tol_U = 0.1;
breaker1 = 0;
breaker2 = 0;
iterations = 0;

Sigma = reshape(experiment.Sigma,[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
Sigma = diag(Sigma);
d = reshape(experiment.d,[],1);

u = u0;
[p,ups] = forward_map(u,params);
data_misfit = norm(Sigma_minus_half * (d - reshape(p,[],1)))^2;
distance_from_u0 = norm(C_minus_half * (u - u0)')^2;
J = data_misfit + distance_from_u0;

tic;

for k = 1:100

    [boldR,curlyR,c] = compute_representers(u,ups,C,params,params.nobtimes);

    h = (u0-u)/(1+alpha) + (boldR * ((curlyR + (1+alpha)*Sigma)\(d - reshape(p,[],1) - c/(1+alpha))))';

    u_cand = u + h;

    [p_cand,ups_cand] = forward_map(u_cand,params);
    data_misfit = norm(Sigma_minus_half * (d - reshape(p_cand,[],1)))^2;
    distance_from_u0 = norm(C_minus_half * (u_cand - u0)')^2;
    J_cand = data_misfit + distance_from_u0;

    iterations = iterations + 1;

    if J_cand < J
        breaker1 = (max(abs(u-u_cand)./abs(u)) < tol_U);
        breaker2 = (J - J_cand < tol_J);
        u = u_cand;
        J = J_cand;
        ups = ups_cand;
        p = p_cand;
        C_k = C - boldR * inv(curlyR + (1+alpha)*Sigma) * boldR';
        
        if plotting
            figure(1)
            disp([J,alpha])
            upper50 = u + 0.674*sqrt(diag(C_k))';
            lower50 = u - 0.674*sqrt(diag(C_k))';
            upper95 = u + 1.96*sqrt(diag(C_k))';
            lower95 = u - 1.96*sqrt(diag(C_k))';
            figure(1)
            fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], 'g', 'EdgeColor','none', 'FaceAlpha',0.25)
            hold on
            plot(params.x_locations,u,"k")
            plot(params.x_locations,upper95,"k--")
            plot(params.x_locations,lower95,"k--")
            plot(experiment.params_fwd.x_locations,experiment.u_true,"r")
            xline(experiment.ups_true(end),"r--")
            hold off
            xlim([0,1])
            ylim([-2,2])
            drawnow
        end
        alpha = alpha/2;
    else
        alpha = alpha*2;
    end

    if breaker1 || breaker2
        u_map = u;
        C_map = C - boldR * inv(curlyR + Sigma) * boldR';
        timer = toc
        break
    end
end

figure(1)
disp([J,alpha])
upper50 = u_map + 0.674*sqrt(diag(C_map))';
lower50 = u_map - 0.674*sqrt(diag(C_map))';
upper95 = u_map + 1.96*sqrt(diag(C_map))';
lower95 = u_map - 1.96*sqrt(diag(C_map))';
figure(1)
fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], 'g', 'EdgeColor','none', 'FaceAlpha',0.25)
hold on
plot(params.x_locations,u_map,"k--")
plot(params.x_locations,upper95,"k")
plot(params.x_locations,lower95,"k")
plot(experiment.params_fwd.x_locations,experiment.u_true,"r")
xline(experiment.ups_true(end),"r--")
hold off
xlim([0,1])
ylim([-2,2])
drawnow

end