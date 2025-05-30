function time_data = find_time_data(obj,t)

times = obj.RTMflow_class.times;
time_data = zeros(t,2);
for time = 1:t
    t_j = obj.physics_class.observation_times(time);
    [~,closest_time] = sort(abs(times-t_j));
    closest_time = closest_time(1);
    time_data(time,1) = min(closest_time+100,length(times)-1);

    del_t = max(times(closest_time+2) - times(closest_time),...
                times(closest_time) - times(closest_time-2));
    time_data(time,2) = del_t^2;
end

end