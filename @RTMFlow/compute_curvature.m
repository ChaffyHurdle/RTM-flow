function curv_edges = compute_curvature(obj,edge_data_t)

    if edge_data_t == 0
        curv_edges = 0;
    else

        front = unique(edge_data_t(:,2:3));
        x = obj.Delaunay_mesh_class.nodes(front,1);
        y = obj.Delaunay_mesh_class.nodes(front,2);
        
        % Step 1: Sort data by y (important for a well-behaved function)
        [y, sortIdx] = sort(y);  
        x = x(sortIdx);
        
        % Step 2: Fit a smoothing spline for x as a function of y
        p = 0.9999; % Smoothing parameter (0 = very smooth, 1 = interpolates all points)
        spline_x = csaps(y, x, p);  
        
        % Step 3: Evaluate spline on a fine grid
        y_fine = linspace(0, 1, 100);
        x_fine = fnval(spline_x, y_fine);
        
        dx_dy = fnval(fnder(spline_x, 1), y_fine); % First derivative
        d2x_dy2 = fnval(fnder(spline_x, 2), y_fine); % Second derivative
        
        % Step 4: Compute curvature
        curvature = 2 * d2x_dy2 ./ (1 + dx_dy.^2).^(3/2); 
        
    
        curv_edges = zeros(1,size(edge_data_t,1));
        for i = 1:size(edge_data_t,1)
            centroid = obj.Delaunay_mesh_class.centroids(edge_data_t(i,1),:);
            [~, dist_id] = min(sum((centroid-[x_fine',y_fine']).^2,2));
            curv_edges(i) = curvature(dist_id);
        end
    end
