function active_sensors = find_active_sensors(obj,t)

times = obj.RTMflow_class.times;
nsensors = obj.physics_class.nsensors;
active_sensors = zeros(t*nsensors,8);
for k = 1:t*nsensors
    % Sensor and time index
    i = obj.i_vec(k);
    j = obj.j_vec(k);

    % Actual time and sensor element index
    t_j = obj.physics_class.observation_times(j);
    x_i = obj.physics_class.sensor_locs(i,:);
    x_i_sensor_elem = obj.RTMflow_class.sensor_element_inds(i);

    % Nearest time to t_j
    time_index = find(obj.RTMflow_class.times > t_j,1);

    % Save data
    is_active = ismember(x_i_sensor_elem, ...
        find(obj.RTMflow_class.all_active_elements(:,time_index)));
    if is_active
        active_sensors(k,1:2) = x_i;
        active_sensors(k,3) = x_i_sensor_elem;
        active_sensors(k,4) = t_j;
        [~,closest_time] = sort(abs(times-t_j));
        closest_time = closest_time(1);
        
        del_t = max(times(closest_time+2) - times(closest_time),...
                    times(closest_time) - times(closest_time-2));

        active_sensors(k,5) = 1;
        active_sensors(k,6) = min(closest_time+100,length(times)-1);
        active_sensors(k,7) = del_t^2;
        active_sensors(k,8) = k;
    end
end

active_sensors = active_sensors(find(active_sensors(:,5)),:);
end