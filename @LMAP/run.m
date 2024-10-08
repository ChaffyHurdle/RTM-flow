function obj = run(obj)

% Converge condition depends on exit criterion for LM algorithm
converged = 0;
iterate = 1;
u = obj.inverse_class.u0;
obj.alpha = 1e5;
obj.tol1 = 0.01;
obj.tol2 = 0.01;

% FWD simulation u0
obj.physics_class.permeability = exp(u);
obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
obj.RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
obj.RTMflow_class = obj.RTMflow_class.run();

% Evaluate posterior estimate at u0
J = obj.evaluate_cost_function(u);

while ~converged && iterate < obj.max_iterations
    disp(iterate)

    % FWD simulation
    obj.physics_class.permeability = exp(u);
    obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
    obj.RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
    obj.RTMflow_class = obj.RTMflow_class.run();

    % Compute lambdas and gradient lambdas
    obj = obj.compute_lambdas();

    % Compute representers
    obj = obj.compute_representers();

    % Compute linearised states
    obj = obj.compute_linearised_states();

    % Construct linear system
    obj = obj.construct_linear_system();

    % Update u_{k} -> u_{k+1}
    obj = obj.compute_h();
    candidate_u = u + h;

    % Check if cost function improved
    obj.physics_class.permeability = exp(candidate_u);
    obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
    obj.RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
    obj.RTMflow_class = obj.RTMflow_class.run();
    candidate_J = obj.evaluate_cost_function(candidate_u);

    % If candidate improves, perform iteration. If not, increase regularisation
    if candidate_J < J
        u = candidate_u;
        J = candidate_J;
        obj.alpha = obj.alpha/10;
    else
        obj.alpha = 10*obj.alpha;
    end
    
    % Convergence check
    if abs(candidate_J - J) < obj.tol1 || sum((u - candidate_u).*obj.mesh_class.element_areas) < obj.tol2
        converged = 1;
    end
    
    iterate = iterate + 1;
end