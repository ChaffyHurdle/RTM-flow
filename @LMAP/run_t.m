function obj = run_t(obj,t)

% Various shorthands
converged = 0;
iterate = 1;
u = obj.u0;
physics_class = obj.physics_class;
mesh_class = obj.mesh_class;
C = obj.inverse_class.C0_inv;
Sigma = diag(reshape(obj.inverse_class.Sigma(:,1:t),[],1));
n_fwds = 5;

% Storage vectors
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
start_time = tic;

% Start loop
while ~converged & iterate < obj.max_iterations

    % Solve adjoint equations, compute representers
    obj = obj.parallel_computations_t(t); % Computes \lambda, \mathbb{R}, \mathcal{R}, d

    % Update u_{k} -> u_{k+1}
    h = obj.compute_h_t(t,n_fwds);
    h_accepted = zeros(1,n_fwds);
    RTM_candidates = cell(1,n_fwds);
    pressure_candidates = cell(1,n_fwds);
    physics_candidates = cell(1,n_fwds);
    J_candidates = zeros(1,n_fwds);
    scaled_misfit_candidates = zeros(1,n_fwds);

    parfor i = 1:n_fwds
        candidate_u = u + h(:,i)';
    
        % Run simulation and check if cost function improved
        candidate_physics = physics_class;
        candidate_physics.permeability = exp(candidate_u');
        candidate_pressure = Pressure(mesh_class,candidate_physics);
        candidate_RTM = RTMFlow(mesh_class,candidate_physics,candidate_pressure,1);
        candidate_RTM = candidate_RTM.run(physics_class.observation_times(t));
        [candidate_J,scaled_data_misfit] = obj.evaluate_cost_function_t(candidate_u,candidate_RTM,t);

        physics_candidates{i} = candidate_physics;
        pressure_candidates{i} = candidate_pressure;
        RTM_candidates{i} = candidate_RTM;
        J_candidates(i) = candidate_J;
        scaled_misfit_candidates(i) = scaled_data_misfit;

        h_accepted(i) = (candidate_J < J);
    end

    if sum(h_accepted) == 0
        time_elapsed = toc(start_time);
        execution_times(iterate+1) = time_elapsed;
        disp("Converged through patience.")
        break
    end

    % Choose step
    firstaccepted = find(h_accepted, 1, 'first');
    %[~,firstaccepted] = min(J_candidates);
    h = h(:,firstaccepted)';
    candidate_u = u + h;
    candidate_RTM = RTM_candidates{firstaccepted};
    candidate_pressure = pressure_candidates{firstaccepted};
    candidate_physics = physics_candidates{firstaccepted};
    candidate_J = J_candidates(firstaccepted);
    scaled_data_misfit = scaled_misfit_candidates(firstaccepted);
    obj.alpha = obj.alpha*obj.scale^(firstaccepted-1);

    disp("Iterate: " + num2str(iterate) + ", Best J: " + num2str(J) + ", Current J: " + num2str(candidate_J) + ", J change: " + num2str(round(100*((candidate_J-J)/J),1)) + "%, Alpha: " + num2str(obj.alpha))
    
    if (abs(candidate_J - J)/abs(J) <= obj.tol2 || max((u - candidate_u)./max(u)) <= obj.tol1) && iterate > 5
        disp("Converged through stopping criterion.")
        converged = 1;
    end

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
    obj.physics_class = candidate_physics;
    obj.pressure_class = candidate_pressure;
    obj.RTMflow_class = candidate_RTM;
    C_post = C - obj.R*((obj.tildePmat + (1+obj.alpha)*Sigma)\obj.R');

    obj.alpha = obj.alpha/obj.scale;
    iterate = iterate + 1;

    % Plot iterate
    obj.plotter(u,C_post,h);
        
end

% Save data 
obj.u_map = u;
obj.C_map = C - obj.R*((obj.tildePmat + Sigma)\obj.R');
obj.u_iterations = u_iterations(:,1:iterate);
obj.J_iterations = J_iterations(1:iterate);
obj.scaled_data_misfit = data_misfit_iterations(1:iterate);
obj.execution_times = [0,execution_times(execution_times>0)];
obj.best_alpha = best_alpha;