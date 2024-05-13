 function opt = cvfem_setup(mesh,bndry,darcy,K)

opt.mesh = mesh;
opt.bndry = bndry;
opt.darcy = darcy;
inlet_location = str2func(opt.bndry.inlet_location_fname);
vent_location = str2func(opt.bndry.vent_location_fname);

NT=size(mesh.elem,1);
c = zeros(NT,2); % this stores a vector with the centers of elements.
for i = 1 : NT
    c(i,:) = mean(mesh.node(mesh.elem(i,:),:));
end

opt.mesh.elem_center=c;
opt.mesh.nnode = size(opt.mesh.node,1);
opt.mesh.nelem = size(opt.mesh.elem,1);

opt = fem2d_init_tri(opt);

opt.mesh.normal_vec = compute_normals(opt.mesh.elem,opt.mesh.node);

opt.cvfem.K = K;

% compute volumes of the control volumes
opt.cvfem.V = compute_volumes(opt.mesh.node,opt.mesh.elem,opt.darcy.thickness);          

[opt.bndry.inlet_flag, opt.bndry.inlet_pos, opt.bndry.Dirichlet] ...
    = inlet_location(opt.mesh.node, opt.mesh.bndry_nodes);
[opt.bndry.vent_flag] ...
    = vent_location(opt.mesh.node,opt.mesh.bndry_nodes);
opt.bndry.vent_idx = find(opt.bndry.vent_flag);

% Nuemann boundary condition
opt.bndry.neumann_flag = find_nuemann_points(opt.mesh.bndry_nodes,opt.bndry.inlet_flag, opt.bndry.vent_flag);


 end