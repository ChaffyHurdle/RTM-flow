function curv = compute_curvature(obj,x)
    
    % Number of points
    n = length(x);
    curv = zeros(n,1); % Initialize curvature array
    
    % Loop over interior points (excluding first and last)
    for i = 2:n-1
        % Get three consecutive points
        x1 = x(i-1,1); y1 = x(i-1,2);
        x2 = x(i,1);   y2 = x(i,2);
        x3 = x(i+1,1); y3 = x(i+1,2);

        % Compute side lengths
        a = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        b = sqrt((x3 - x2)^2 + (y3 - y2)^2);
        c = sqrt((x3 - x1)^2 + (y3 - y1)^2);
        
        % Compute signed triangle area
        A = 0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)); 

        % Compute circumradius (avoid division by zero)
        R = (a * b * c) / (4 * abs(A)); % Absolute area for radius
        cross_product = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2);
        curv(i) = 2 / R * sign(cross_product); % Assign sign based on cross product

    end
end
