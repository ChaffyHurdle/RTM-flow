function  [vent_flag, vent_elem] = vent_location(node,bnd_node)
% Function identifies nodes that lie on the inlet boundary and currently
% returns a zero for elements that lie on the inlet boundary???

% Sets the outlet at y==0

% number of mesh nodes
nnode = size(node,1);

%number of boundary nodes
nbnd_node = nnz(bnd_node); 

vent = zeros(nbnd_node,1);
vent_idx = 0;

candidate = find(bnd_node);

for i = 1 : nbnd_node

    ci = candidate(i);

    if node(ci,2) == 0.1  
        vent_idx = vent_idx+1;
        vent(vent_idx) = ci;
    end
end
vent_flag = sparse(nnode,1);
vent_flag(vent(1:vent_idx)) = 1;

vent_elem = 0;

end