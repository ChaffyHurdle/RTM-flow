function obj = plot_sensitivity(obj)

nodes = obj.mesh_class.nodes;
elements = obj.mesh_class.elements;
num_ob_times = length(obj.darcy_class_u.observation_times);
edge_data = obj.all_edge_data;
is_moving_boundary_ob_times = obj.is_moving_boundary_ob_times;
is_moving_boundary_ob_times_u_plus_h = obj.is_moving_boundary_ob_times_u_plus_h;

figure
for i = 1:num_ob_times
    free_nodes = obj.active_nodes_u(:,i) & ~obj.is_moving_boundary_ob_times(:,i);
    max_diff = max(max(abs(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*free_nodes),...
        max(abs(obj.pressures_u_plus_h(:,i) - obj.pressures_u(:,i) - obj.p_tildes(:,i)).*free_nodes));

    % |p_{u+h} - p_{u}|
    subplot(2,num_ob_times,i)
    pdeplot(nodes',elements', ...
        XYData=abs(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*free_nodes, ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([0,max_diff])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$|p_{u+h} - p_u|$", 'interpreter', 'latex')

    % |p_{u+h} - (p_u + \tilde{p})|
    subplot(2,num_ob_times,i+num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=abs(obj.pressures_u_plus_h(:,i) - obj.pressures_u(:,i) - obj.p_tildes(:,i)).*free_nodes, ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([0,max_diff])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$|p_{u+h} - (p_u + \tilde{p})|$", 'interpreter', 'latex')
end

figure
for i = 1:num_ob_times
    subplot(4,num_ob_times,i)
    free_nodes = obj.active_nodes_u(:,i) & ~obj.is_moving_boundary_ob_times(:,i);
    max_cbar = max(max(obj.p_tildes(:,i).*free_nodes),max( (obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*free_nodes) );
    min_cbar = min(min(obj.p_tildes(:,i).*free_nodes),min( (obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*free_nodes) );

    pdeplot(nodes',elements', ...
        XYData=(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*free_nodes, ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$p_{u+h}-p_u$",'interpreter','latex')

    subplot(4,num_ob_times,i+num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=obj.p_tildes(:,i).*free_nodes, ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$\tilde{p}$",'interpreter','latex')

    subplot(4,num_ob_times,i+2*num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=(obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i)).*is_moving_boundary_ob_times(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$p_{u+h}-p_u$ on boundary",'interpreter','latex')

    subplot(4,num_ob_times,i+3*num_ob_times)
    pdeplot(nodes',elements', ...
        XYData=obj.p_tildes(:,i).*is_moving_boundary_ob_times(:,i), ...
        XYStyle='flat',ColorMap="jet",Mesh="off")
    caxis([min_cbar,max_cbar])
    hold on
    plot(nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),1),nodes(boolean(is_moving_boundary_ob_times_u_plus_h(:,i)),2),'w.')
    hold off
    hold on
    for j = 1:length(edge_data{obj.time_inds_u(i)})
        plot([nodes(edge_data{obj.time_inds_u(i)}(j,2),1), nodes(edge_data{obj.time_inds_u(i)}(j,3),1)],...
            [nodes(edge_data{obj.time_inds_u(i)}(j,2),2), nodes(edge_data{obj.time_inds_u(i)}(j,3),2)],'w');
    end
    hold off
    title("$\tilde{p}$ on boundary",'interpreter','latex')
end

figure
for i = 1:num_ob_times
    plot(i, sum( (obj.pressures_u_plus_h(:,i)-obj.pressures_u(:,i) - obj.p_tildes(:,i)).^2 )/sum(free_nodes),'k*' )
    hold on
end
hold off

obj.p_tildes_at_sensors = zeros(size(obj.RTMflow_u.pressure_data));
for i = 1:width(obj.p_tildes)
    p_tilde = obj.p_tildes(:,i);
    p_tilde_at_sensors = interpolate_pressures(p_tilde,obj.RTMflow_u,obj.mesh_class);
    obj.p_tildes_at_sensors(:,i) = p_tilde_at_sensors;
end


end

%% The following functions are designed to interpolate pressure at sensor locations
function pressure_at_sensors = interpolate_pressures(p_tilde,obj,mesh)

    sensor_element_inds = obj.sensor_element_inds;
    sensor_locs = obj.physics_class.sensor_locs;
    pressure_at_sensors = zeros(1,length(sensor_locs));

    for i = 1:length(sensor_locs)
        sensor_loc = sensor_locs(i,:);
        nodes_surrounding_inds = mesh.elements(sensor_element_inds(i),:);
        nodes_surrounding = mesh.nodes(nodes_surrounding_inds,:);
        [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, sensor_loc);
        P1 = p_tilde(nodes_surrounding_inds(1));
        P2 = p_tilde(nodes_surrounding_inds(2));
        P3 = p_tilde(nodes_surrounding_inds(3));
        pressure_at_sensors(i) = lambda1 * P1 + lambda2 * P2 + lambda3 * P3;
    end
end

% Auxiliary function (previously defined)
function [lambda1, lambda2, lambda3] = compute_barycentric_coords(nodes_surrounding, p_tilde)
    v1 = nodes_surrounding(1,:); v2 = nodes_surrounding(2,:); v3 = nodes_surrounding(3,:);
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));
    lambda1 = ((v2(2) - v3(2))*(p_tilde(1) - v3(1)) + (v3(1) - v2(1))*(p_tilde(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p_tilde(1) - v3(1)) + (v1(1) - v3(1))*(p_tilde(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end