function num_taken = cpu_scheduler(obj_data,active_sensors)

% Number of cores/workers
pool = gcp();
numWorkers = pool.NumWorkers;

% Sort the active inds
active_inds = find(active_sensors);
active_times = obj_data.j_vec(active_inds);
[~,sorted_active_times_inds] = sort(active_times,'descend');
sorted_inds = active_inds(sorted_active_times_inds);

temp_sorted_inds = active_times(sorted_active_times_inds);
% for i =1:5
%     disp(sum(temp_sorted_inds==i))
% end
num_taken = 0;
total_taken = 0;
count = 0;
while total_taken ~= length(sorted_inds)
    count = count + 1;
    if isscalar(unique(temp_sorted_inds)) || length(temp_sorted_inds) <= numWorkers
        num_taken = [num_taken length(temp_sorted_inds)];
        total_taken = total_taken+length(temp_sorted_inds);
    else
        first_element = temp_sorted_inds(1);
        num_first_element = sum(temp_sorted_inds == first_element);
        num_to_take = numWorkers*ceil(num_first_element/numWorkers);
        num_taken = [num_taken num_to_take];
        total_taken = total_taken + num_to_take;
        temp_sorted_inds = temp_sorted_inds(num_to_take+1:end);
    end
    if count > 100
        disp("Count loop forever")
        break
    end
end
% Fast correction for an open debugging issue
final_ind = find(num_taken, 1, 'last');
num_taken = num_taken(1:final_ind);