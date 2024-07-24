function obj = run(obj)

% Converge condition depends on exit criterion for LM algorithm
converged = 0;
iterate = 1;
u = obj.inverse_class.u0;

while ~converged
    disp(iterate)

    % FWD simulation
    obj.physics_class.permeability = exp(u);
    obj.pressure_class = Pressure(obj.mesh_class,obj.physics_class);
    obj.RTMflow_class = RTMFlow(obj.mesh_class,obj.physics_class,obj.pressure_class);
    obj.RTMflow_class.visualise_class.is_plotting_volume = true;
    obj.RTMflow_class = obj.RTMflow_class.run();

    % Compute lambdas
    obj = obj.compute_lambdas();

    % Compute representers
    obj = obj.compute_representers();

    % Compute linearised states

    % Construct linear system

    % Update u_{k} -> u_{k+1}
    u = u + normrnd(0,5,1,length(u));

    iterate = iterate + 1;
    if iterate == 2
        converged = 1;
    end
end