function obj = solve(obj)

while ~obj.is_fully_saturated()

    %% Solve Pressure Problem
    obj.pressure_class = obj.pressure_class.solve();
    
    %% Solve flow problem
    obj = obj.compute_flow_rates();
    
    %% Increment to new time
    obj = obj.compute_time_increment();
    obj.time = obj.time + obj.time_step;

    %% Update flow volumes and moving boundaries
    obj = obj.update_filling_percentage();
    obj = obj.update_computational_domain();

end

disp("end")

end