function [U,timer,total_iterations] = run_RML(u0,C,params,experiment,plotting)

samples = 5000;
prior_samples = mvnrnd(u0,C,samples);
U = zeros(params.Nx,samples);
total_iterations = 0;

C_minus_half = inv(sqrtm(C));
alpha = 1e4;
tol_J = 0.1;
tol_U = 0.1;
breaker1 = 0;
breaker2 = 0;

Sigma = reshape(experiment.Sigma,[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
Sigma = diag(Sigma);
original_d = reshape(experiment.d,[],1);

delete(gcp('nocreate'))
i=10;
fprintf('Number of slots available: %d\n', i);
parpool('Threads', i);

tic;

parfor i = 1:samples
    u0 = prior_samples(i,:);
    u = u0;
    d = original_d + normrnd(reshape(experiment.Sigma,[],1)*0,reshape(experiment.Sigma,[],1));
    [p,ups] = forward_map(u,params);
    data_misfit = norm(Sigma_minus_half * (d - reshape(p,[],1)))^2;
    distance_from_ui = norm(C_minus_half * (u - u0)')^2;
    J = data_misfit + distance_from_ui;
    alpha = 1e4;
    tol_J = 0.1;
    tol_U = 0.1;
    breaker1 = 0;
    breaker2 = 0;
    iteration_count = 0;

    for j = 1:100
        [boldR,curlyR,c] = compute_representers(u,ups,C,params,5);
        h = (u0-u)/(1+alpha) + (boldR * ((curlyR + (1+alpha)*Sigma)\(d - reshape(p,[],1) - c/(1+alpha))))';
    
        u_cand = u + h;
    
        [p_cand,ups_cand] = forward_map(u_cand,params);
        data_misfit = norm(Sigma_minus_half * (d - reshape(p_cand,[],1)))^2;
        distance_from_ui = norm(C_minus_half * (u_cand - u0)')^2;
        J_cand = data_misfit + distance_from_ui;
        iteration_count = iteration_count+1;

        if J_cand < J
            breaker1 = (max(abs(u-u_cand)./abs(u)) < tol_U);
            breaker2 = (J - J_cand < tol_J);
            u = u_cand;
            J = J_cand;
            ups = ups_cand;
            p = p_cand;
            alpha = alpha/2;
        else
            alpha = alpha*2;
        end

        if breaker1 || breaker2
            U(:,i) = u;
            total_iterations = total_iterations + iteration_count;
            break
        end
    end
    
end
timer = toc

if plotting
    figure(3)
    U_mean=mean(U,2);
    X=[params.x_locations,fliplr(params.x_locations)];
    Y=[(U_mean-0.674*sqrt(var(U,0,2)))',fliplr((U_mean+0.674*sqrt(var(U,0,2)))')];
    fill(X,Y, 'g', 'EdgeColor','none', 'FaceAlpha',0.25);
    hold on
    plot(params.x_locations,U_mean,'k--')
    plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
    plot(params.x_locations,U_mean+1.96*sqrt(var(U,0,2)),'k')
    plot(params.x_locations,U_mean-1.96*sqrt(var(U,0,2)),'k')
    xline(experiment.ups_true(end),"r--")
    %legend('$$U_{EKI}\pm 0.674\sigma_{EKI}$$','$$U_{EKI}$$','truth','location','north','fontsize',20,'interpreter','latex')
    xlim([0,1])
    ylim([-2,2])
    hold off
    drawnow
end

end