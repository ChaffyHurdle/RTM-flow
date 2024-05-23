function [cvfem, bndry] = update_comp_domain(cvfem,mesh,bndry)

%% unpacking
activeElement = cvfem.activeElement;
activeNode = cvfem.activeNode;
elem = mesh.elem;
nnode = mesh.nnode;
new_filled_volume = cvfem.new_filled_volume;
has_node_i = mesh.has_node_i;
neumann_flag = bndry.neumann_flag;
fFactor = cvfem.fFactor;

nelem = size(elem,1);
newActiveElement = zeros(nelem,1);
candidate_elem = zeros(nelem,1);
bnd_node = zeros(nnode,1);
Dirichlet = zeros(nnode,1);

% Update the activeElement array
for i = 1 : length(new_filled_volume)
    ni = new_filled_volume(i);
    for j = 2 : has_node_i(ni,1)+1
        candidate_elem(has_node_i(ni,j)) = 1;
    end
end
newActiveElement(activeElement<1 & candidate_elem ==1) = 1;
elems = elem(newActiveElement==1,:);
activeNode(elems(:))=1;
activeElement = min(ones(nelem,1),activeElement + newActiveElement);

% aelem_idx = find(activeElement==1);
aelem = elem(activeElement==1,:);

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
        if sum(activeElement(has_node_i(bnd_idx_i,2:has_node_i(bnd_idx_i,1)+1)))<has_node_i(bnd_idx_i,1)
            Dirichlet(bnd_idx_i) = 1;
        end
    end
end
Dirichlet(fFactor>0&fFactor<1) = 1;

%% Repacking
cvfem.activeElement = activeElement;
cvfem.activeNode = activeNode;
cvfem.newActiveElement = newActiveElement;
bndry.Dirichlet = Dirichlet;

end
