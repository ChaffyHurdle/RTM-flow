function [cvfem , bndry]= still_solver(cvfem,mesh,bndry)

%% build/update FEM stiffness matrix
cvfem.A = assembling_stiffness_tri(mesh.node, mesh.elem, cvfem.K, cvfem.A, cvfem.newActiveElement);

%% impose dirichlet boundary conditions
[u, b, freeNode, bdNode] = bcond_dirichlet(cvfem.A,bndry.Dirichlet,cvfem.activeNode,mesh.node,bndry);
b = modify_rhs(cvfem.A,u,b);

%% solve matrix equation for pressure
cvfem.u = solve(cvfem.A,u,b,freeNode);
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