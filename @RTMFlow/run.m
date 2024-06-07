function obj = run(obj)

while ~obj.is_fully_saturated()

    %% Solve pressure & velocity Problem
    obj.pressure_class = obj.pressure_class.solve();
    obj.velocity_class = obj.velocity_class.compute_velocity(obj.pressure_class);
    
    %% Solve flow problem
    obj = obj.compute_flow_rates();

    %% Visualise
    obj.visualise_class.plot(obj);
    obj.visualise_class.animate(obj);
    
    %% Increment to new time
    obj = obj.update_time_level();
 
    %% Update flow volumes and moving boundaries
    obj = obj.update_filling_percentage();
    obj = obj.update_computational_domain();

    %% Compute any voids/vacuum
    obj = obj.find_dryspots();
    obj = obj.void_partition();
    obj = obj.apply_ideal_gas_law();

end

disp("end")

end