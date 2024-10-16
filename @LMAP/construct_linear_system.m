function obj = construct_linear_system(obj)

nobs = length(obj.physics_class.observation_times);
nsensors = obj.physics_class.nsensors;

tildePmat = zeros(nobs*nsensors,nobs*nsensors);
%tildePmat2 = zeros(nobs*nsensors,nobs*nsensors);

% ob_time_inds = [];
% for i = 1:nobs
%     t_ind = find(obj.RTMflow_class.times > obj.physics_class.observation_times(i),1)-1;
%     ob_time_inds = [ob_time_inds t_ind];
% end

% for k = 1:nobs*nsensors
%     p_tilde_i = obj.p_tildes(:,ob_time_inds,k);
%     p_tilde_i_sensors = zeros(nsensors,nobs);
%     for i = 1:nobs
%         p_tilde_i_sensors(:,i) = interpolate_ptildes(p_tilde_i(:,i),obj);
%     end
%     tildePmat(:,k) = reshape(p_tilde_i_sensors',1,[]);
% end

% for j = 1:nobs*nsensors
%     p_tilde_j = obj.p_tildes(:,:,j);
%     for k = 1:nobs*nsensors
%         tildePmat(k,j) = interpolate_ptildes(p_tilde_j,k,obj);
%     end
% end

for j = 1:nobs*nsensors
    for k = j:nobs*nsensors
        tildePmat(j,k) = sum(obj.Q(:,j).*obj.R(:,k).*obj.mesh_class.element_areas);
        tildePmat(k,j) = tildePmat(j,k);
    end
end

d = zeros(nsensors*nobs,1);
for k = 1:nobs*nsensors
    %d(k) = interpolate_ptildes(obj.p_tilde0,k,obj);
    d(k) = sum(obj.Q(:,k).*(obj.u0 - obj.u)'.*obj.mesh_class.element_areas);
end
obj.tildePmat = tildePmat;
%obj.tildePmat2 = tildePmat2;
obj.d = d;
end


%% The following functions are designed to interpolate pressure at sensor locations
function ptilde_at_sensor = interpolate_ptildes(ptilde,k,obj)

    i = obj.i_vec(k);
    j = obj.j_vec(k);

    time_index = find(obj.RTMflow_class.times > obj.physics_class.observation_times(j),1)-1;
    ptilde_t = ptilde(:,time_index);
    sensor_element = obj.RTMflow_class.sensor_element_inds(i);
    sensor_loc = obj.physics_class.sensor_locs(i,:);

    mesh = obj.mesh_class;

    nodes_surrounding_inds = mesh.elements(sensor_element,:);
    nodes_surrounding = mesh.nodes(nodes_surrounding_inds,:);
    [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, sensor_loc);
    P1 = ptilde_t(nodes_surrounding_inds(1));
    P2 = ptilde_t(nodes_surrounding_inds(2));
    P3 = ptilde_t(nodes_surrounding_inds(3));
    ptilde_at_sensor = lambda1 * P1 + lambda2 * P2 + lambda3 * P3;

end

% Auxiliary function (previously defined)
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p(1) - v3(1)) + (v3(1) - v2(1))*(p(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p(1) - v3(1)) + (v1(1) - v3(1))*(p(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end
