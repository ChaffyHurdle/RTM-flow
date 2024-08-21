function obj = find_sensor_locs_mesh(obj)

sensor_locs = obj.physics_class.sensor_locs;

% Find nearest sensor locations on computational mesh (indices)
obj.sensor_inds_on_mesh = dsearchn(obj.Delaunay_mesh_class.nodes,sensor_locs);

% Locate sensors on computational mesh
obj.sensor_locs_on_mesh = obj.Delaunay_mesh_class.nodes(obj.sensor_inds_on_mesh,:);

sensor_element_inds = zeros(1,length(sensor_locs));

for i = 1:length(sensor_locs)
    sensor_loc = sensor_locs(i,:);
    sensor_element_inds(i) = findElementContainingPoint(obj.Delaunay_mesh_class,sensor_loc);
end
obj.sensor_element_inds = sensor_element_inds;


function elementIndex = findElementContainingPoint(mesh, point)
    
    nodes = mesh.nodes;
    elements = mesh.elements;
    elementIndex = NaN;

    for j = 1:size(elements, 1)

        elementNodes = elements(j, :);
        v1 = nodes(elementNodes(1), :);
        v2 = nodes(elementNodes(2), :);
        v3 = nodes(elementNodes(3), :);

        % Compute the barycentric coordinates for the point
        [lambda1, lambda2, lambda3] = computeBarycentricCoords(v1, v2, v3, point);

        % Check if the point is inside the triangle (all barycentric coordinates between 0 and 1)
        if lambda1 >= 0 && lambda2 >= 0 && lambda3 >= 0 && ...
           lambda1 <= 1 && lambda2 <= 1 && lambda3 <= 1
            elementIndex = j;
            return;
        end
    end
end

function [lambda1, lambda2, lambda3] = computeBarycentricCoords(v1, v2, v3, p)
    % Compute the denominator of the barycentric coordinates
    denominator = (v2(2) - v3(2))*(v1(1) - v3(1)) + (v3(1) - v2(1))*(v1(2) - v3(2));

    % Compute the barycentric coordinates
    lambda1 = ((v2(2) - v3(2))*(p(1) - v3(1)) + (v3(1) - v2(1))*(p(2) - v3(2))) / denominator;
    lambda2 = ((v3(2) - v1(2))*(p(1) - v3(1)) + (v1(1) - v3(1))*(p(2) - v3(2))) / denominator;
    lambda3 = 1 - lambda1 - lambda2;
end


end