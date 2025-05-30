function load_vector = compute_load_vec(obj,active_sensors_t,t_ind)

nodes = obj.mesh_class.nodes;
elements = obj.mesh_class.elements;
times = obj.RTMflow_class.times;

load_vector = zeros(obj.mesh_class.num_nodes,size(active_sensors_t,1));
for k = 1: size(active_sensors_t,1)
    x_i = active_sensors_t(k,1:2);
    x_i_sensor_elem = active_sensors_t(k,3);
    t_j = active_sensors_t(k,4);
    dt = active_sensors_t(k,7);

    % Compute aspect of load vector from \delta(x_i-x) contribution.
    nodes_surrounding_xi = nodes(elements(x_i_sensor_elem,:),:);
    [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding_xi, x_i);
    load_vector(elements(x_i_sensor_elem,1),k) = lambda1;
    load_vector(elements(x_i_sensor_elem,2),k) = lambda2;
    load_vector(elements(x_i_sensor_elem,3),k) = lambda3;
    load_vector(:,k) = load_vector(:,k)*obj.delta_t(t_j, ...
                                                    times(t_ind), ...
                                                    dt);
end

end

%% Compute barycentric coords
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p(1) - v3(1)) + (v3(1) - v2(1))*(p(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p(1) - v3(1)) + (v1(1) - v3(1))*(p(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end