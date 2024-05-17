%% Problem set up
my_mesh = Mesh(p,e,t);
my_darcy = Darcy(0.1, 0.35, 0.2, @permeability);
my_pressure = Pressure(my_mesh,@is_inlet,@is_vent,@p_D);
my_volume = Volume(my_pressure,my_darcy);
%my_visuals = Visualisation();

%% compile CVFEM method
my_cvfem = CVFEM(my_mesh,my_pressure,my_volume,my_darcy,[],[]);

%% Execute solver
my_cvfem.solve()

%% Argument set up
function K = permeability(x)
    
    K = [1e-10 0; 0 1e-10];
end

function p = p_D(points,inlet_nodes,outlet_nodes,time)

p = zeros(size(points,1),1);
p(inlet_nodes) = 1.5e5;

end

function bool = is_inlet(node)

bool = (node(2) == 0);

end

function bool = is_vent(node)

bool = (node(2) == 1);

end
    
function [inlet_flag, inlet_pos, Dirichlet] = inlet_location(node,bnd_node)

nnode = size(node,1);
nbnd_node = nnz(bnd_node);

inlet = zeros(nbnd_node,1);
inlet_idx = 0;

candidate = find(bnd_node);
for i = 1 : nbnd_node
    
    ci = candidate(i);
    if  node(ci,2) == 0
        inlet_idx = inlet_idx+1;
        inlet(inlet_idx) = ci;
    end
end

inlet_flag = sparse(nnode,1);
inlet_flag(inlet(1:inlet_idx)) = 1;

inlet_pos = node(inlet_flag==1,:);

% sets Dirichlet condition
Dirichlet =zeros(nnode,1);
Dirichlet(inlet_flag==1) = 1;

end

function  [vent_flag, vent_elem] = vent_location(node,bnd_node)

nnode = size(node,1);
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





    