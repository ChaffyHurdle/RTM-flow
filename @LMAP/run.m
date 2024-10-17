function obj = run(obj)

% Initialise quantities
converged = 0;
iterate = 1;
do_over = 0;
u = obj.u0;
patience = 0;
C = obj.inverse_class.C0_inv;
Sigma = diag(reshape(obj.inverse_class.Sigma,[],1));
G_u0 = obj.RTMflow_class.pressure_data;
u_iterations = zeros(obj.mesh_class.num_elements,obj.max_iterations);
J_iterations = zeros(1,obj.max_iterations);
execution_times = zeros(1,obj.max_iterations);
best_alpha = obj.alpha;

% Evaluate posterior cost function at u0
J = obj.evaluate_cost_function(u,obj.RTMflow_class);
u_iterations(:,1) = u;
J_iterations(1) = J;
execution_times(1) = 0;


while ~converged & iterate < obj.max_iterations

    if patience > 3
        disp("Converged through patience")
        break
    end

    % Compute lambdas, gradient lambdas, representers and linearised states
    if ~do_over
        start_time = tic;
        obj = obj.parallel_computations();
        obj = obj.construct_linear_system(); % compute R, d
    end

    % Update u_{k} -> u_{k+1}
    h = obj.compute_h();
    candidate_u = u + h;

    % Run simulation and check if cost function improved
    candidate_physics = obj.physics_class;
    candidate_physics.permeability = exp(candidate_u');
    candidate_pressure = Pressure(obj.mesh_class,candidate_physics);
    candidate_RTM = RTMFlow(obj.mesh_class,candidate_physics,candidate_pressure);
    candidate_RTM = candidate_RTM.run();
    candidate_J = obj.evaluate_cost_function(candidate_u,candidate_RTM);

    disp("Iterate: " + num2str(iterate) + ", Best J: " + num2str(J) + ", Current J: " + num2str(candidate_J) + ", J change: " + num2str(round(100*((candidate_J-J)/J),1)) + "%, Alpha: " + num2str(obj.alpha))

    % If candidate improved, perform iteration. If not, increase regularisation
    if candidate_J < J
        % Convergence check
        if abs(candidate_J - J)/abs(J) <= obj.tol1 && iterate>3 %|| max((u - candidate_u)')/max(u') < obj.tol2
            disp("LMAP converged")
            converged = 1;
        end

        disp("Cost function improved")
        % Save data
        u_iterations(:,iterate+1) = candidate_u;
        J_iterations(iterate+1) = candidate_J;
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
            plot(obj.inverse_class.data(i,:),'k*')
            hold on
            plot(candidate_RTM.pressure_data(i,:),'r*')
            plot(G_u0(i,:),'b*')
            hold off
            ylim([0.9,2.1])
            title("$D$ (black), G($u_0$) (blue), G($u_{map}$) (red)",'interpreter','latex')
        end
        drawnow

        posterior_samples = mvnrnd(u,C_post,25);
        figure(3)
        for i = 1:25
            subplot(5,5,i)
            pdeplot(obj.mesh_class.nodes',...
                    obj.mesh_class.elements', ...
                    XYData=posterior_samples(i,:), ...
                    XYStyle='interp',ColorMap="jet",Mesh="off")
            clim([c_min,c_max])
        end
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
obj.execution_times = execution_times(1:iterate);
obj.best_alpha = best_alpha;