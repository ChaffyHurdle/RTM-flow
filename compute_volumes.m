function V  = compute_volumes(node,elem,thickness)
N = size(node,1);
NT = size(elem,1);

V = zeros(N,1);
for i = 1 : NT
    elem_i = elem(i,1:3);
    V(elem_i) = V(elem_i) + compute_volume_i(node(elem_i,:), thickness);
end

end

%% computes volume on an individual element
function v = compute_volume_i(xi, thickness)

A = .5*det([1 xi(1,1) xi(1,2);...
         1 xi(2,1) xi(2,2);...
         1 xi(3,1) xi(3,2)]);
         
v = ones(3,1)*A/3*thickness;


end


