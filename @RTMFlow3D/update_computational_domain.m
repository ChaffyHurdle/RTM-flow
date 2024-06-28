function obj = update_computational_domain(obj)

%% unpacking
active_elements = obj.active_elements;
is_node_active = obj.pressure_class.is_node_active;
new_filled_volume = obj.new_filled_volume;
Dirichlet = 0*obj.pressure_class.is_Dirichlet;
fFactor = obj.volume_fill_percentage;

%% add new element to the volume
active_elements(new_filled_volume) = 1;
newActiveElement = false(obj.Delaunay_mesh_class.num_elements,1);
newActiveElement(new_filled_volume) = true;
new_active_nodes = obj.Delaunay_mesh_class.elements(newActiveElement,:);
new_active_nodes = new_active_nodes(:);
is_node_active(new_active_nodes) = true;

%% update Dirichlet boundary conditions on moving front
active_element_faces = obj.Delaunay_mesh_class.element_faces(active_elements>=0.5,:);

%% determine uniqueness of active faces
sorted_faces = sort(active_element_faces,2);
[u,~,n] = unique(sorted_faces,'rows');
counts = accumarray(n(:), 1);
sorted_exteriorF = u(counts == 1,:);
boundary_faces = ismember(sorted_faces,sorted_exteriorF,'rows');
boundary_nodes = active_element_faces(boundary_faces,:);
boundary_nodes = unique(boundary_nodes(:));
          
for i = 1:length(boundary_nodes)

    node_num = boundary_nodes(i);

    if ~obj.pressure_class.is_Neumann(node_num)

        Dirichlet(node_num)=1;
    end

end


elem = obj.Delaunay_mesh_class.elements;
a = elem(fFactor>0,:); b = elem(fFactor<1,:);
Dirichlet(unique(intersect(a(:),b(:)))) = 1;


%{
aelem = elem(active_elements==1,:);

total_edge = [aelem(:,[2 3]); aelem(:,[3,1]); aelem(:,[1,2])];


Cnode = sparse(total_edge(:,1), total_edge(:,2),1,nnode,nnode);
[r,c] = find(Cnode-Cnode');
bnd_node([r;c]) = 1;

bnd_idx = find(bnd_node>0);

for i = 1:length(bnd_idx)
    bnd_idx_i = bnd_idx(i);
    if ~neumann_flag(bnd_idx_i)
        Dirichlet(bnd_idx_i)=1;
    else % We make the intersection between flow front and wall Dirichelt
        if sum(active_elements(has_node_i(bnd_idx_i,2:has_node_i(bnd_idx_i,1)+1)))<has_node_i(bnd_idx_i,1)
            Dirichlet(bnd_idx_i) = 1;
        end
    end
end
Dirichlet(fFactor>0&fFactor<1) = 1;
%}

%% Repacking
obj.pressure_class.is_node_active = is_node_active;
obj.active_elements = active_elements;
obj.pressure_class.new_active_elements = newActiveElement;
obj.pressure_class.is_Dirichlet = Dirichlet;

%% set Dirichlet BC values
obj.pressure_class.pressure = obj.pressure_class.p_D(obj.pressure_class);

end