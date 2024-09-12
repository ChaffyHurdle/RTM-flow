%% Clearing and loading
clear;
close all;
clc;
delete(gcp('nocreate'))
addpath('Meshes')
meshes = {'p.mat', 'e.mat', 't.mat'};
for i = 1:numel(meshes)
    load(meshes{i})
end

%% Mesh set up
my_mesh = DelaunayMesh(p,e,t);

var_matern = 0.25; length_scale = 0.1; nu_matern = 1.5;
matern_args = [var_matern,length_scale,nu_matern];
my_inverse = Inversion(my_mesh,my_mesh,matern_args);

%% Generate u and h
my_inverse = my_inverse.generate_u();
my_inverse.plot_u_true();
h = my_inverse.generate_h();

%% Darcy setup
mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;

% Approx. ob. times for 5 equal increments of the front (hard-coded to work 
% for mean 0 prior generating u_true). Stop when ~95% filled.
observation_times = linspace(0.1,0.9,5).^2*mu*phi/(2*(p_I-p_0));
T = 0.95^2*mu*phi/(2*(p_I-p_0));

% Set sensor locs
sensor_locs_x = [0.25,0.5,0.75];
sensor_locs_y = [0.25,0.5,0.75];
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

% Define true permeability and place within physics class
K_true = exp(my_inverse.u_true);
my_darcy = Physics(mu, phi, thickness, p_I, p_0, K_true, sensor_locs, observation_times,T);

%% Sensitivity setup
my_sensitivity = Sensitivity(my_mesh,my_darcy,my_inverse.u_true,h/10);
my_sensitivity = my_sensitivity.fwd_solves();
my_sensitivity = my_sensitivity.run();

for i = 1:5
    [ d, ix_u_plus_h ] = min( abs( my_sensitivity.RTMflow_u_plus_h.times-my_darcy.observation_times(i) ) );
    [ d, ix_u ] = min( abs( my_sensitivity.RTMflow_u.times-my_darcy.observation_times(i) ) );
    max_diff = max(max(abs(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h)-my_sensitivity.RTMflow_u.pressures(:,ix_u))),...
        max(abs(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h) - my_sensitivity.RTMflow_u.pressures(:,ix_u) - my_sensitivity.p_tilde(:,ix_u))));
    figure(i)

    % p_{u+h}
    subplot(2,3,1)
    pdeplot(my_mesh.nodes',...
                my_mesh.elements', ...
                XYData=my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    title("$p_{u+h}$", 'interpreter', 'latex')
    % p_{u}
    subplot(2,3,2)
    pdeplot(my_mesh.nodes',...
                my_mesh.elements', ...
                XYData=my_sensitivity.RTMflow_u.pressures(:,ix_u).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    title("$p_u$", 'interpreter', 'latex')
    % |p_{u+h} - p_{u}|
    subplot(2,3,3)
    pdeplot(my_mesh.nodes',...
                my_mesh.elements', ...
                XYData=abs(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h)-my_sensitivity.RTMflow_u.pressures(:,ix_u)).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    caxis([0,max_diff])
    title("$|p_{u+h} - p_u|$", 'interpreter', 'latex')
    % p_{u+h}
    subplot(2,3,4)
    pdeplot(my_mesh.nodes',...
                my_mesh.elements', ...
                XYData=my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    title("$p_{u+h}$", 'interpreter', 'latex')
    % p_{u} + \tilde{p}
    subplot(2,3,5)
    pdeplot(my_mesh.nodes',...
                my_mesh.elements', ...
                XYData=(my_sensitivity.RTMflow_u.pressures(:,ix_u) + my_sensitivity.p_tilde(:,ix_u)).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    title("$p_u + \tilde{p}$", 'interpreter', 'latex')
    % |p_{u+h} - (p_u + \tilde{p})|
    subplot(2,3,6)
    pdeplot(my_mesh.nodes',...
            my_mesh.elements', ...
            XYData=abs(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h) - my_sensitivity.RTMflow_u.pressures(:,ix_u) - my_sensitivity.p_tilde(:,ix_u)).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    caxis([0,max_diff])
    title("$|p_{u+h} - (p_u + \tilde{p})|$", 'interpreter', 'latex')
    ob_time = "t_" + num2str(i);
    sgtitle(ob_time)
end

for i = 1:5
    [ d, ix_u_plus_h ] = min( abs( my_sensitivity.RTMflow_u_plus_h.times-my_darcy.observation_times(i) ) );
    [ d, ix_u ] = min( abs( my_sensitivity.RTMflow_u.times-my_darcy.observation_times(i) ) );
    max_cbar = max(max(my_sensitivity.p_tilde(:,ix_u)),max(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h)-my_sensitivity.RTMflow_u.pressures(:,ix_u)));
    min_cbar = min(min(my_sensitivity.p_tilde(:,ix_u)),min(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h)-my_sensitivity.RTMflow_u.pressures(:,ix_u)));
    figure(i+5)
    subplot(1,2,1)
    pdeplot(my_mesh.nodes',...
                    my_mesh.elements', ...
                    XYData=(my_sensitivity.RTMflow_u_plus_h.pressures(:,ix_u_plus_h)-my_sensitivity.RTMflow_u.pressures(:,ix_u)).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    caxis([min_cbar,max_cbar])
    title("$p_{u+h}-p_u$",'interpreter','latex')
    subplot(1,2,2)
    pdeplot(my_mesh.nodes',...
                    my_mesh.elements', ...
                    XYData=my_sensitivity.p_tilde(:,ix_u).*my_sensitivity.RTMflow_u.active_nodes(:,ix_u),XYStyle='flat',ColorMap="jet",Mesh="on")
    caxis([min_cbar,max_cbar])
    title("$\tilde{p}$",'interpreter','latex')
    ob_time = "t_" + num2str(i);
    sgtitle(ob_time)
end

%% Plot differences in moving front and pressure

figure(1)
k = 1;
for i = my_darcy.observation_times
    subplot(1,5,k)
    u_RTMflow = my_sensitivity.RTMflow_u;
    [ d, ix ] = min( abs( u_RTMflow.times-i ) );
    front = boolean(u_RTMflow.active_nodes(:,ix) & u_RTMflow.Dirichlet_nodes(:,ix) & ~u_RTMflow.pressure_class.is_inlet);
    front_nodes = sortrows(my_mesh.nodes(front,:),2);
    plot(front_nodes(:,1),front_nodes(:,2))
    
    hold on
    u_plus_h_RTMflow = my_sensitivity.RTMflow_u_plus_h;
    [ d, ix ] = min( abs( u_plus_h_RTMflow.times-i ) );
    front = boolean(u_plus_h_RTMflow.active_nodes(:,ix) & u_plus_h_RTMflow.Dirichlet_nodes(:,ix) & ~u_plus_h_RTMflow.pressure_class.is_inlet);
    front_nodes = sortrows(my_mesh.nodes(front,:),2);
    plot(front_nodes(:,1),front_nodes(:,2))
    hold off

    xlim([0,1]);ylim([0,1])
    k = k+1;
end

figure(2)
u_RTMflow = my_sensitivity.RTMflow_u;
max_error = zeros(1,4);
k = 1;
for i=1:4
    u_plus_h_RTMflow = my_sensitivity.RTMflow_u_plus_h;
    for j=my_darcy.observation_times
        subplot(4,5,k)
        [ d, ix_u_plus_h ] = min( abs( u_plus_h_RTMflow.times-j ) );
        [ d, ix_u ] = min( abs( u_RTMflow.times-j ) );
        h_norm = max(2^(-i)*my_sensitivity.h);
        pdeplot(my_mesh.nodes',...
               my_mesh.elements', ...
               XYData=(u_plus_h_RTMflow.pressures(:,ix_u_plus_h) - u_RTMflow.pressures(:,ix_u))/h_norm,XYStyle='interp',ColorMap="jet",Mesh="off")
        colorbar; caxis([-1/2 1/2]);
        k = k + 1;
    end
    max_error(i) = (u_plus_h_RTMflow.pressures(500,ix_u_plus_h) - u_RTMflow.pressures(500,ix_u))/h_norm;
end

for i = 1:50
    figure(1)
    subplot(1,2,1)
    pdeplot(my_mesh.nodes',...
            my_mesh.elements', ...
            XYData=u_RTMflow.all_active_elements(:,i),XYStyle='flat',ColorMap="jet",Mesh="on")
    title(num2str(sum(u_RTMflow.active_nodes(:,i))))
    subplot(1,2,2)
    pdeplot(my_mesh.nodes',...
            my_mesh.elements', ...
            XYData=u_RTMflow.all_active_elements(:,i+1),XYStyle='flat',ColorMap="jet",Mesh="on")
    title(num2str(sum(u_RTMflow.active_nodes(:,i+1))))
     input("Another? Y/Q");
end

moving_bndry = double(u_RTMflow.active_nodes & u_RTMflow.Dirichlet_nodes & ~u_RTMflow.pressure_class.is_inlet);
for i = 1:50
    figure(1)
    pdeplot(my_mesh.nodes',...
            my_mesh.elements', ...
            XYData=moving_bndry(:,i),XYStyle='flat',ColorMap="jet",Mesh="on")
end

%% Edges
t = 25;
active_at_t = my_sensitivity.RTMflow_u.all_active_elements(:,t);
active_nodes = boolean(my_sensitivity.RTMflow_u.active_nodes(:,t)); % nodes in D(t)
dirichlet_nodes = boolean(my_sensitivity.RTMflow_u.Dirichlet_nodes(:,t)); % inlet, front and beyond front
inlet_nodes = my_sensitivity.RTMflow_u.pressure_class.is_inlet;
is_moving_boundary = dirichlet_nodes ...
                   & active_nodes ...
                   & ~inlet_nodes;
moving_boundary_inds = find(is_moving_boundary);

boundaryEdges = [];
boundaryTriangles = [];
for i = 1:length(moving_boundary_inds)
    
    nodei = moving_boundary_inds(i);

    % Step 3: Find all triangles that include this boundary node
    [rowIdx, ~] = find(my_mesh.elements == nodei);  % Find row indices where the node is used in T
    
    % Step 4: Loop over each triangle containing this boundary node
    for j = 1:length(rowIdx)
        triangleIdx = rowIdx(j);  % Get the triangle index
        nodes = my_mesh.elements(triangleIdx, :); % Get the nodes of the triangle
        
        % Step 5: Check if the edges of the triangle are on the boundary
        n1 = nodes(1); n2 = nodes(2); n3 = nodes(3);

        if sum(is_moving_boundary(nodes)) < 3
            % Check edge {n1, n2}
            if is_moving_boundary(n1) == 1 && is_moving_boundary(n2) == 1 && active_nodes(n3)
                boundaryEdges = [boundaryEdges; n1, n2];
                boundaryTriangles = [boundaryTriangles; triangleIdx];
            end
            
            % Check edge {n2, n3}
            if is_moving_boundary(n2) == 1 && is_moving_boundary(n3) == 1 && active_nodes(n1)
                boundaryEdges = [boundaryEdges; n2, n3];
                boundaryTriangles = [boundaryTriangles; triangleIdx];
            end
            
            % Check edge {n3, n1}
            if is_moving_boundary(n3) == 1 && is_moving_boundary(n1) == 1 && active_nodes(n2)
                boundaryEdges = [boundaryEdges; n3, n1];
                boundaryTriangles = [boundaryTriangles; triangleIdx];
            end
        
        end
        
    end
end
boundaryEdges = unique(sort(boundaryEdges, 2), 'rows');
boundaryTriangles = unique(boundaryTriangles);

figure(1)
subplot(1,2,1)
pdeplot(my_mesh.nodes',...
        my_mesh.elements', ...
        XYData=active_at_t,XYStyle='flat',ColorMap="jet",Mesh="on")
subplot(1,2,2)
triangle_loc = zeros(length(my_mesh.elements),1);
triangle_loc(boundaryTriangles) = 1;
pdeplot(my_mesh.nodes',...
        my_mesh.elements', ...
        XYData=triangle_loc,XYStyle='flat',ColorMap="jet",Mesh="on")
hold on
for i = 1:length(boundaryEdges)
    plot([my_mesh.nodes(boundaryEdges(i,1),1),my_mesh.nodes(boundaryEdges(i,2),1)], ...
        [my_mesh.nodes(boundaryEdges(i,1),2),my_mesh.nodes(boundaryEdges(i,2),2)],'w')
end
hold off

%% CV linearised states
RTMflow_u = my_sensitivity.RTMflow_u;
tilde_p = zeros(size(RTMflow_u.pressures));
v_h = zeros(size(RTMflow_u.pressures));

for t = 200:300
    boundary_nodes_t = RTMflow_u.active_nodes(:,t) & RTMflow_u.Dirichlet_nodes(:,t);
    moving_boundary_nodes_t = boundary_nodes_t & ~RTMflow_u.pressure_class.is_inlet;

    boundary_nodes_t_plus = RTMflow_u.active_nodes(:,t+1) & RTMflow_u.Dirichlet_nodes(:,t+1);
    moving_boundary_nodes_t_plus = boundary_nodes_t_plus & ~RTMflow_u.pressure_class.is_inlet;

    new_node = find(abs(boundary_nodes_t_plus-boundary_nodes_t));
    
    boundary_node_inds = find(boundary_nodes_t);
    candidate_elem = zeros(my_mesh.num_elements,1);

    for i = 1 : length(boundary_node_inds)
        candidate = boundary_node_inds(i);
        for j = 2 : RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,1)+1
            elem = RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,j);
            candidate_elem(RTMflow_u.Voronoi_mesh_class.has_node_i(candidate,j)) = 1;
        end
    end
    % candidate_elem = find(candidate_elem);

    boundaryEdges = [];
    boundaryTriangles = [];
    for i = 1:length(boundary_node_inds)
        
        nodei = boundary_node_inds(i);
    
        % Step 3: Find all triangles that include this boundary node
        [rowIdx, ~] = find(my_mesh.elements == nodei);  % Find row indices where the node is used in T
        
        % Step 4: Loop over each triangle containing this boundary node
        for j = 1:length(rowIdx)
            triangleIdx = rowIdx(j);  % Get the triangle index
            nodes = my_mesh.elements(triangleIdx, :); % Get the nodes of the triangle
            
            % Step 5: Check if the edges of the triangle are on the boundary
            n1 = nodes(1); n2 = nodes(2); n3 = nodes(3);
    
            if sum(boundary_nodes_t(nodes)) < 3
                % Check edge {n1, n2}
                if boundary_nodes_t(n1) == 1 && boundary_nodes_t(n2) == 1 && RTMflow_u.active_nodes(n3,t)
                    boundaryEdges = [boundaryEdges; n1, n2];
                    boundaryTriangles = [boundaryTriangles; triangleIdx];
                end
                
                % Check edge {n2, n3}
                if boundary_nodes_t(n2) == 1 && boundary_nodes_t(n3) == 1 && RTMflow_u.active_nodes(n1,t)
                    boundaryEdges = [boundaryEdges; n2, n3];
                    boundaryTriangles = [boundaryTriangles; triangleIdx];
                end
                
                % Check edge {n3, n1}
                if boundary_nodes_t(n3) == 1 && boundary_nodes_t(n1) == 1 && RTMflow_u.active_nodes(n2,t)
                    boundaryEdges = [boundaryEdges; n3, n1];
                    boundaryTriangles = [boundaryTriangles; triangleIdx];
                end
            
            end
            
        end
    end
    boundaryEdges = unique(sort(boundaryEdges, 2), 'rows');

    figure(1)
    pdeplot(my_mesh.nodes',...
            my_mesh.elements', ...
            XYData=candidate_elem,XYStyle='flat',ColorMap="jet",Mesh="on")
    hold on
    scatter(my_mesh.nodes(moving_boundary_nodes_t,1),my_mesh.nodes(moving_boundary_nodes_t,2),'wo')
    for i = 1:length(boundaryEdges)
        plot([my_mesh.nodes(boundaryEdges(i,1),1),my_mesh.nodes(boundaryEdges(i,2),1)], ...
            [my_mesh.nodes(boundaryEdges(i,1),2),my_mesh.nodes(boundaryEdges(i,2),2)],'w')
    end
    hold off
    input('')
end 
