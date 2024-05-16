function Q = update_flow_rate_tri(opt)

Q = zeros(opt.mesh.nnode,1);
candidate_elem = zeros(opt.mesh.nelem,1);
u = opt.cvfem.u;
K = opt.cvfem.K;
darcy = opt.darcy;
fFactor = opt.cvfem.fFactor;
normal_vec = opt.mesh.normal_vec;

%% Allocating velocity field
velocity = zeros(opt.mesh.nelem,2);

candidate_node = find(~xor(opt.cvfem.activeNode==1,1-fFactor > 1e-15));

for i = 1 : length(candidate_node)
    candidate = candidate_node(i);
    for j = 2 : opt.mesh.has_node_i(candidate,1)+1
        elem = opt.mesh.has_node_i(candidate,j);
        if opt.cvfem.activeElement(elem) >= 0.5
            candidate_elem(opt.mesh.has_node_i(candidate,j)) = 1;
        end
    end
end

candidate_elem = find(candidate_elem);

%% partitions elements into inlet-lying and non inlet lying elements 
% H.W: could be simplified here
inlet_elem = candidate_elem(opt.mesh.elem_including_inlet_edge(candidate_elem)==1);
other_elem = candidate_elem(opt.mesh.elem_including_inlet_edge(candidate_elem)==0);

%% Computes volume fluxes across elements connected to the inlet boundary
for i = 1 : length(inlet_elem)
    elem_i = opt.mesh.elem(inlet_elem(i),:);
    node_i = opt.mesh.node(elem_i,:);
    vi = velocity_centre_tri(u(elem_i),node_i,K,darcy, darcy.phi(inlet_elem(i)));
    velocity(inlet_elem(i),:) = vi;
    Q(elem_i) = Q(elem_i) + local_flux_tri_inlet(node_i,vi,fFactor(elem_i),opt.bndry.inlet_flag(elem_i),opt.mesh.bndry_nodes(elem_i),darcy);
end

%% Computes volume fluxes across all remaining elements
for i = 1 : length(other_elem)
    elem_idx = other_elem(i);
    elem_i = opt.mesh.elem(elem_idx,:);
    vi = velocity_centre_tri(u(elem_i),opt.mesh.node(elem_i,:),K,darcy, darcy.phi(elem_idx));
    velocity(elem_idx,:) = vi;
    Q(elem_i) = Q(elem_i) + local_flux_tri(normal_vec(elem_idx,:),vi,darcy);
end
end

function velocity = velocity_centre_tri(p,tri_nodes,K,darcy_para, phi)

%% Element geometric properties
x = tri_nodes(:,1); y = tri_nodes(:,2);
centroid = mean(tri_nodes);
area = 0.5*abs(x'*circshift(y,-1) - circshift(x,-1)'*y);

%% Evaluate pressure gradient
grad_phi = [y(2)-y(3) y(3)-y(1) y(1)-y(2); x(3)-x(2) x(1)-x(3) x(2)-x(1)]/2/area;
pressure_gradient = grad_phi * p;

%% Apply Darcy's Law (-K/(mu*phi) * grad (p)
visocity = darcy_para.mu;
porosity = phi;
velocity = -K(centroid)*pressure_gradient/(visocity*porosity);
end

%% function to compute outflow of standard elements
function qn = local_flux_tri(normal_vec,v,darcy_para)
qn = zeros(3,1);

qn(1) = normal_vec(:,[1 2])*v;
qn(2) = normal_vec(:,[3 4])*v;
qn(3) = normal_vec(:,[5 6])*v;

qn = qn*darcy_para.thickness;

end

%% function to compute inflow/outflow of boundary-lying elements
function qn = local_flux_tri_inlet(x,v,fFactor,inlet_flag, bnd_flag,darcy_para)
if sum(bnd_flag)==3
    error('This code does not support an triangular element with three boundary nodes.');
end
qn = zeros(3,1);

o = sum(x)/3;
a = sum(x([1 2],:))/2;
b = sum(x([2 3],:))/2;
c = sum(x([1 3],:))/2;

n1 = [o(2) - a(2), a(1)-o(1)];
n2 = [o(2) - b(2), b(1)-o(1)];
n3 = [o(2) - c(2), c(1)-o(1)];


if sum(inlet_flag) == 2 % The element contains an edge lying on the inlet
    if inlet_flag(1) == 0
        n4 = fliplr(.5*x(2,:)-.5*x(3,:)).*[1 -1];
        temp1 = zeros(1,2);
        temp2 = n4-n2;
        temp3 = n4+n2;
        if sum(1-fFactor < eps)>= 1
            temp1 = temp1-n1+n3;
            temp2 = temp2 + n1;
            temp3 = temp3 - n3;
        end
    elseif inlet_flag(2) == 0
        n4 = fliplr(.5*x(3,:)-.5*x(1,:)).*[1 -1];
        temp1 = n4+n3;
        temp2 = zeros(1,2);
        temp3 = n4-n3;
        if sum(1-fFactor < eps)>= 1
            temp1 = temp1-n1;
            temp2 = temp2 + n1 - n2;
            temp3 = temp3 - n2;
        end
    elseif inlet_flag(3) == 0
        n4 = fliplr(.5*x(1,:)-.5*x(2,:)).*[1 -1];
        temp1 = n4-n1;
        temp2 = n4+n1;
        temp3 = zeros(1,2);
        if sum(1-fFactor < eps)>= 1
            temp1 = temp1+n3;
            temp2 = temp2-n2;
            temp3 = temp3 +n2- n3;
        end
    end
    
    qn(1) = temp1*v;
    qn(2) = temp2*v;
    qn(3) = temp3*v;
      
elseif sum(inlet_flag) == 1 % At this moment, we ignore elements having single inlet node.
    if sum(1-fFactor < eps)>= 1
        qn(1) = -(n1-n3)*v;
        qn(2) = -(n2-n1)*v;
        qn(3) = -(n3-n2)*v;
    end
else
    error('The triangular element must have at most two nodes on the inlet boundary.');
end

qn = qn*darcy_para.thickness;

end