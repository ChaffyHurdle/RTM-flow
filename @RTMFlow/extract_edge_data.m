function edgedata = extract_edge_data(obj)

active_nodes = obj.pressure_class.is_node_active; % nodes in D(t)
Dirichlet_nodes = obj.pressure_class.is_Dirichlet;
is_moving_boundary = active_nodes & Dirichlet_nodes & ~obj.pressure_class.is_inlet;
moving_boundary_inds = find(is_moving_boundary);
mesh_class = obj.Delaunay_mesh_class;

% Find elements that make up the moving boundary and their connectivity
edgedata = [];
for i = 1:length(moving_boundary_inds)
    
    nodei = moving_boundary_inds(i);
    [rowIdx, ~] = find(mesh_class.elements == nodei);
    
    for j = 1:length(rowIdx)
        triangleIdx = rowIdx(j);
        nodes = mesh_class.elements(triangleIdx, :);
       
        n1 = nodes(1); n2 = nodes(2); n3 = nodes(3);

        if sum(is_moving_boundary(nodes)) < 3
            if is_moving_boundary(n1) == 1 && is_moving_boundary(n2) == 1 && active_nodes(n3)
                edgedata = [edgedata; triangleIdx n1 n2];
            end
            
            if is_moving_boundary(n2) == 1 && is_moving_boundary(n3) == 1 && active_nodes(n1)
                edgedata = [edgedata; triangleIdx n2 n3];
            end
            
            if is_moving_boundary(n3) == 1 && is_moving_boundary(n1) == 1 && active_nodes(n2)
                edgedata = [edgedata; triangleIdx n3 n1];
            end
        
        end
        
    end
end
% Remove duplicate rows 
edgedata = unique(edgedata,'rows');

% Compute unit outer normals 
edges = mesh_class.nodes(edgedata(:,2),:) - mesh_class.nodes(edgedata(:,3),:);
normals = [-edges(:,2), edges(:,1)];
normals = normals./vecnorm(normals')';
neg_rows = normals(:, 1) < 0;
normals(neg_rows, :) = -normals(neg_rows, :);
edgedata = [edgedata normals];

end