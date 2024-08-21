function obj = run(obj)

mesh = obj.Delaunay_mesh_class;
observation_times = obj.physics_class.observation_times;
sensor_inds_mesh = obj.sensor_inds_on_mesh;
t_index = 1;
p_old = obj.pressure_class.pressure;
p_old_at_sensors = interpolate_pressures(p_old,obj,mesh);
p_gradients = obj.pressure_class.pressure_gradient;
Q_old = obj.volume_rates_of_flow;
filling_facs_old = obj.volume_fill_percentage;
stiffness = obj.pressure_class.stiffness_matrix;
Dirichlets = obj.pressure_class.is_Dirichlet;
active_nodes = obj.pressure_class.is_node_active;
active_elements = obj.active_elements;
new_active_elements = obj.pressure_class.new_active_elements;
t_old = obj.time;
it = 1;

obj = obj.add_data_all_times(it,t_old,p_old,...
    p_gradients,Q_old,filling_facs_old,stiffness,active_nodes,Dirichlets, ...
    active_elements,new_active_elements);

tic

while ~obj.is_fully_saturated() && obj.time + obj.time_step <= obj.physics_class.T

    it = it + 1;

    %% Solve pressure & velocity problem
    obj.pressure_class = obj.pressure_class.solve();
    obj.velocity_class = obj.velocity_class.compute_velocity(obj.pressure_class);
    
    %% Solve flow problem
    obj = obj.compute_flow_rates();

    %% Visualise
    % obj.visualise_class.plot(obj);
    
    %% Increment to new time
    obj = obj.update_time_level();
    t_new = obj.time;
    dt = obj.time_step;
 
    %% Update flow volumes and moving boundaries
    obj = obj.update_filling_percentage();

    %% Save data to object (each time)
    p_new = obj.pressure_class.pressure;
    p_new_at_sensors = interpolate_pressures(p_new,obj,mesh);
    Q_new = obj.volume_rates_of_flow;
    filling_facs_new = obj.volume_fill_percentage;
    p_gradients_new = obj.pressure_class.pressure_gradient;
    stiffness = obj.pressure_class.stiffness_matrix;
    Dirichlets = obj.pressure_class.is_Dirichlet;
    active_nodes = obj.pressure_class.is_node_active;
    active_elements = obj.active_elements;
    new_active_elements = obj.pressure_class.new_active_elements;
    obj = obj.add_data_all_times(it,t_new,p_new,...
        p_gradients_new,Q_new,filling_facs_new,...
        stiffness,active_nodes,Dirichlets,active_elements,new_active_elements);

    %% Update domain
    obj = obj.update_computational_domain();

    
    %% Save pressure at sensors
    if t_index <= length(observation_times)
        if (observation_times(t_index) > t_old) && (observation_times(t_index) < t_new)
            observation_time = observation_times(t_index);
            disp([t_old,observation_time,t_new])
            obj.pressure_data(:,t_index) ...
                = p_old_at_sensors + ...
                ((observation_time - t_old)/dt)*(p_new_at_sensors - p_old_at_sensors);
            t_index = t_index + 1;
            obj.visualise_class.plot(obj);
        end

    %else
    %    break
    end
    
    p_old = p_new;
    p_old_at_sensors = p_new_at_sensors;
    t_old = t_new;
end

obj.wall_time = toc;
disp("Wall-time elapsed: " + num2str(obj.wall_time) + ' s')

if t_index <= length(observation_times)
    obj.pressure_data(:,t_index:end) = p_new(sensor_inds_mesh);
end

end

%% The following functions are designed to interpolate pressure at sensor locations
function pressure_at_sensors = interpolate_pressures(p,obj,mesh)

    sensor_element_inds = obj.sensor_element_inds;
    sensor_locs = obj.physics_class.sensor_locs;
    pressure_at_sensors = zeros(1,length(sensor_locs));

    for i = 1:length(sensor_locs)
        sensor_loc = sensor_locs(i,:);
        nodes_surrounding_inds = mesh.elements(sensor_element_inds(i),:);
        nodes_surrounding = mesh.nodes(nodes_surrounding_inds,:);
        [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, sensor_loc);
        P1 = p(nodes_surrounding_inds(1));
        P2 = p(nodes_surrounding_inds(2));
        P3 = p(nodes_surrounding_inds(3));
        pressure_at_sensors(i) = lambda1 * P1 + lambda2 * P2 + lambda3 * P3;
    end
end

% Auxiliary function (previously defined)
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p(1) - v3(1)) + (v3(1) - v2(1))*(p(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p(1) - v3(1)) + (v1(1) - v3(1))*(p(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end
