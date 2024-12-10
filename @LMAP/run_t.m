function obj = run_t(obj,t)

% Initialise quantities
converged = 0;
iterate = 1;
do_over = 0;
u = obj.u0;
patience = 0;

C = obj.inverse_class.C0_inv;
Sigma = diag(reshape(obj.inverse_class.Sigma(:,1:t),[],1));
Sigma_minus_half = diag(1./sqrt(reshape(obj.inverse_class.Sigma(:,1:t),[],1)));
G_u0 = obj.RTMflow_class.pressure_data(:,1:t);
u_iterations = zeros(obj.mesh_class.num_elements,obj.max_iterations);
J_iterations = zeros(1,obj.max_iterations);
data_misfit_iterations = zeros(1,obj.max_iterations);
execution_times = zeros(1,obj.max_iterations);
best_alpha = obj.alpha;

% Evaluate posterior cost function at u0
[J,scaled_misfit_J] = obj.evaluate_cost_function_t(u,obj.RTMflow_class,t);
u_iterations(:,1) = u;
J_iterations(1) = J;
data_misfit_iterations(1) = scaled_misfit_J;
execution_times(1) = 0;


while ~converged & iterate < obj.max_iterations

    if patience > 3
        disp("Converged through patience")
        break
    end

    dscrpncy = sum((Sigma_minus_half * (reshape(obj.inverse_class.data(:,1:t),[],1) - reshape(obj.RTMflow_class.pressure_data(:,1:t),[],1))).^2);
    disp([dscrpncy,chi2inv(0.95,numel(obj.inverse_class.Sigma(:,1:t)))])
    if dscrpncy < chi2inv(0.95,numel(obj.inverse_class.Sigma(:,1:t)))
        disp("Converged through discrepancy")
        break
    end

    % Compute lambdas, gradient lambdas, representers and linearised states
    if ~do_over
        start_time = tic;
        obj = obj.parallel_computations_t(t);
        obj = obj.construct_linear_system_t(t); % compute R, d
    end

    % Update u_{k} -> u_{k+1}
    h = obj.compute_h_t(t);
    candidate_u = u + h;

    % Run simulation and check if cost function improved
    candidate_physics = obj.physics_class;
    candidate_physics.permeability = exp(candidate_u');
    candidate_pressure = Pressure(obj.mesh_class,candidate_physics);
    candidate_RTM = RTMFlow(obj.mesh_class,candidate_physics,candidate_pressure);
    candidate_RTM = candidate_RTM.run(inf);
    [candidate_J,scaled_data_misfit] = obj.evaluate_cost_function_t(candidate_u,candidate_RTM,t);

    disp("Iterate: " + num2str(iterate) + ", Best J: " + num2str(J) + ", Current J: " + num2str(candidate_J) + ", J change: " + num2str(round(100*((candidate_J-J)/J),1)) + "%, Alpha: " + num2str(obj.alpha))

    % If candidate improved, perform iteration. If not, increase regularisation
    if candidate_J < J
        % Convergence check
        if (abs(candidate_J - J)/abs(J) <= obj.tol1 || max((u - candidate_u))/max(u) <= obj.tol2) && iterate > 5 %|| max((u - candidate_u)')/max(u') < obj.tol2
            disp("LMAP converged")
            converged = 1;
        end

        disp("Cost function improved")
        % Save data
        u_iterations(:,iterate+1) = candidate_u;
        J_iterations(iterate+1) = candidate_J;
        data_misfit_iterations(iterate+1) = scaled_data_misfit;
        time_elapsed = toc(start_time);
        execution_times(iterate+1) = time_elapsed;

        best_alpha = min(best_alpha,obj.alpha);

        % Various sets
        u = candidate_u;
        obj.u = u;
        J = candidate_J;
        obj.alpha = obj.alpha/2;
        obj.physics_class = candidate_physics;
        obj.pressure_class = candidate_pressure;
        obj.RTMflow_class = candidate_RTM;
        C_post = C - obj.R*inv(obj.tildePmat + (1+obj.alpha)*Sigma)*obj.R';
        do_over = 0;
        patience = 0;
        iterate = iterate + 1;

        % Plot iterate
        figure(1)
        subplot(2,2,1)
        c_min = min(min(obj.inverse_class.u_true),min(u));
        c_max = max(max(obj.inverse_class.u_true),max(u));
        pdeplot(obj.inverse_class.fwd_mesh.nodes',...
                obj.inverse_class.fwd_mesh.elements', ...
                XYData=obj.inverse_class.u_true, ...
                XYStyle='interp',ColorMap="jet",Mesh="off")
        clim([c_min,c_max])
        title("True u")
        hold on
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
        hold off
        subplot(2,2,2)
        pdeplot(obj.mesh_class.nodes',...
                obj.mesh_class.elements', ...
                XYData=u, ...
                XYStyle='interp',ColorMap="jet",Mesh="off")
        clim([c_min,c_max])
        title("$u_{map}$",'interpreter','latex')
        hold on
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
        hold off
        subplot(2,2,3)
        pdeplot(obj.mesh_class.nodes',...
                obj.mesh_class.elements', ...
                XYData=h, ...
                XYStyle='interp',ColorMap="jet",Mesh="off")
        hold on
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
        hold off
        title("h")
        subplot(2,2,4)
        pdeplot(obj.mesh_class.nodes',...
                obj.mesh_class.elements', ...
                XYData=diag(C_post), ...
                XYStyle='interp',ColorMap="jet",Mesh="off")
        clim([0,0.25])
        hold on
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
        hold off
        title("$C_{map}$",'interpreter','latex')
        drawnow

        figure(2)
        for i = 1:length(obj.physics_class.sensor_locs)
            subplot(sqrt(length(obj.physics_class.sensor_locs)),sqrt(length(obj.physics_class.sensor_locs)),i)
            plot(obj.inverse_class.data(i,:),'ko')
            hold on
            plot(candidate_RTM.pressure_data(i,:),'r*')
            plot(G_u0(i,:),'b*')
            hold off
            ylim([obj.physics_class.p_0-0.1,obj.physics_class.p_I+0.1])
            title("$D$ (black), G($u_0$) (blue), G($u_{map}$) (red)",'interpreter','latex')
        end

        % posterior_samples = mvnrnd(u,C_post,25);
        % figure(3)
        % for i = 1:25
        %     subplot(5,5,i)
        %     pdeplot(obj.mesh_class.nodes',...
        %             obj.mesh_class.elements', ...
        %             XYData=posterior_samples(i,:), ...
        %             XYStyle='interp',ColorMap="jet",Mesh="off")
        %     clim([c_min,c_max])
        % end


        figure(4)
        centroid_x = obj.inverse_class.inv_mesh.centroids(:,1);
        centroid_y = obj.inverse_class.inv_mesh.centroids(:,2);
        
        % Compute distance-based transparency
        % Initialize transparency as ones (fully opaque)
        alpha_map = diag(C_post)/obj.inverse_class.matern_var;
        
        % Create an alpha overlay
        F = scatteredInterpolant(centroid_x, centroid_y, u', 'linear', 'none');
        [xq, yq] = meshgrid(linspace(0,1,1000), ...
                            linspace(0,1,1000));
        zq = F(xq, yq);
        
        % Create an alpha overlay
        F = scatteredInterpolant(centroid_x, centroid_y, alpha_map, 'linear', 'none');
        [xq, yq] = meshgrid(linspace(0,1,1000), ...
                            linspace(0,1,1000));
        alpha_overlay = F(xq, yq);
        
        % Plot the alpha overlay as an image
        im = imagesc(linspace(0,1,1000), linspace(0,1,1000), zq);
        colormap jet
        set(im, 'AlphaData', 1 - alpha_overlay); % Transparency decreases with alpha_map
        set(im, 'AlphaDataMapping', 'none'); % Prevent scaling of alpha
        clim([c_min,c_max])
        colorbar
        set(gca,'YDir','normal')
        xlim([0,1])
        ylim([0,1])
        hold on
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'wo','filled')
        scatter(obj.physics_class.sensor_locs(:,1),obj.physics_class.sensor_locs(:,2),'ko')
        hold off
        drawnow
        
    else
        disp("Cost function got worse, increasing regularisation")
        obj.alpha = 2*obj.alpha;
        do_over = 1;
        patience = patience + 1;
    end
end

obj.u_map = u;
obj.C_map = C - obj.R*inv(obj.tildePmat + Sigma)*obj.R';
obj.u_iterations = u_iterations(:,1:iterate);
obj.J_iterations = J_iterations(1:iterate);
obj.scaled_data_misfit = data_misfit_iterations(1:iterate);
obj.execution_times = execution_times(1:iterate);
obj.best_alpha = best_alpha;