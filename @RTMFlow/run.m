function obj = run(obj)

observation_times = obj.observation_times;
sensor_inds_mesh = obj.sensor_inds_on_mesh;
t_index = 1;
p_old = obj.pressure_class.pressure;
Q_old = obj.volume_rates_of_flow;
filling_facs_old = obj.volume_fill_percentage;
t_old = obj.time;

obj.times = [obj.times t_old];
obj.pressures = [obj.pressures p_old];
obj.pressure_gradients{1} = obj.pressure_class.pressure_gradient;
obj.flow_rates = [obj.flow_rates Q_old];
obj.filling_factors = [obj.filling_factors filling_facs_old];
it = 2;

while ~obj.is_fully_saturated()

    %% Solve pressure & velocity problem
    obj.pressure_class = obj.pressure_class.solve();
    p_new = obj.pressure_class.pressure;
    obj.pressures = [obj.pressures p_new];
    obj.pressure_gradients{it} = obj.pressure_class.pressure_gradient;
    it = it + 1;
    obj.velocity_class = obj.velocity_class.compute_velocity(obj.pressure_class);
    
    %% Solve flow problem
    obj = obj.compute_flow_rates();
    Q_new = obj.volume_rates_of_flow;
    obj.flow_rates = [obj.flow_rates Q_new];

    %% Visualise
    obj.visualise_class.plot(obj);
    
    %% Increment to new time
    obj = obj.update_time_level();
 
    %% Update flow volumes and moving boundaries
    obj = obj.update_filling_percentage();
    filling_facs_new = obj.volume_fill_percentage;
    obj.filling_factors = [obj.filling_factors filling_facs_new];
    obj = obj.update_computational_domain();

    t_new = obj.time;
    obj.times = [obj.times t_new];
    dt = obj.time_step;
    
    % Save pressure at sensors
    if t_index <= length(observation_times)
        if (observation_times(t_index) > t_old) && (observation_times(t_index) < t_new)
            observation_time = observation_times(t_index);
            disp([t_old,observation_time,t_new])
            obj.pressure_data(:,t_index) ...
                = p_old(sensor_inds_mesh) + ((observation_time - t_old)/dt)*(p_new(sensor_inds_mesh) - p_old(sensor_inds_mesh));
            t_index = t_index + 1;
            %obj.visualise_class.plot(obj);
        end

    %else
    %    break
    end
    
    p_old = p_new;
    t_old = t_new;
end

disp("end")
if t_index <= length(observation_times)
    obj.pressure_data(:,t_index:end) = p_new(sensor_inds_mesh);
end

end