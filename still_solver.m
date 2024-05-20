function [cvfem , bndry]= still_solver(cvfem,mesh,bndry)

[cvfem, bndry] = update_comp_domain(cvfem,mesh,bndry);

%% build/update FEM stiffness matrix
cvfem.A = assembling_stiffness_tri(mesh.node, mesh.elem, cvfem.K, cvfem.A, cvfem.newActiveElement);

%% impose dirichlet boundary conditions
[u, b, freeNode, bdNode] = bcond_dirichlet(cvfem.A,bndry.Dirichlet,cvfem.activeNode,mesh.node,bndry);
b = modify_rhs(cvfem.A,u,b);

%% solve matrix equation for pressure
cvfem.u = solve(cvfem.A,u,b,freeNode);
end

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

function [u, b, freeNode, bdNode] = bcond_dirichlet(A,Dirichlet,activeNode, node, bndry)
g_D = str2func(bndry.gd_filename);
N = size(A,1);
b = zeros(N,1);
isBdNode = false(N,1);
isBdNode(Dirichlet==1) = true;
bdNode = find(isBdNode);
freeNode = find(~isBdNode);
freeNode = freeNode(activeNode(freeNode)==1);
u = zeros(N,1);

u(bdNode) = g_D(node(bdNode,:),bndry);
end

function b = modify_rhs(A,u,b)
b = b - A*u;
end

function u = solve(A,u,b,freeNode)
u(freeNode) = A(freeNode,freeNode)\b(freeNode);
end