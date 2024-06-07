function obj = compute_flow_rates(obj)

%% Extract needed information
fFactor = obj.volume_fill_percentage;
normal_vec = obj.volume_class.volume_outflow_vectors;

%% Sets up flow rates if needed
obj.volume_rates_of_flow = zeros(obj.mesh_class.num_nodes,1);

%% construct a list of possible nodes/elements to be added
candidate_elem = zeros(obj.mesh_class.num_elements,1);
candidate_node = find(~xor(obj.is_node_active==1,1-obj.volume_fill_percentage > 1e-15));

for i = 1 : length(candidate_node)
    candidate = candidate_node(i);
    for j = 2 : obj.volume_class.has_node_i(candidate,1)+1
        elem = obj.volume_class.has_node_i(candidate,j);
        if obj.active_elements(elem) >= 0.5
            candidate_elem(obj.volume_class.has_node_i(candidate,j)) = 1;
        end
    end
end
candidate_elem = find(candidate_elem);

for i = 1:length(candidate_elem)

    %% extract local element properties
    element = obj.mesh_class.elements(candidate_elem(i),:);
    nodes = obj.mesh_class.nodes(element,:);

    %% extract velocity in element centre
    vi = obj.velocity_class.velocity(candidate_elem(i),:)';

    is_inlet_connected = sum(obj.inlet_flag(element)) ~= 0;

    if is_inlet_connected
        
        bnd_flag = obj.inlet_flag;
        bnd_flag(obj.mesh_class.boundary_nodes) = 1;
        
        local_flow_rate = local_flux_tri_inlet(nodes,vi,fFactor(element),obj.inlet_flag(element),bnd_flag(element),obj.darcy_class);
    else

        local_flow_rate = local_flux_tri(normal_vec(candidate_elem(i),:),vi,obj.darcy_class);
    end

    obj.volume_rates_of_flow(element) = obj.volume_rates_of_flow(element)...
                                                        + local_flow_rate;

end

end


%% function to compute outflow of standard elements
function qn = local_flux_tri(normal_vec,v,darcy_class)
qn = zeros(3,1);

qn(1) = normal_vec(:,[1 2])*v;
qn(2) = normal_vec(:,[3 4])*v;
qn(3) = normal_vec(:,[5 6])*v;

qn = qn*darcy_class.thickness;

end

function qn = local_flux_tri_inlet(x,v,fFactor,inlet_flag, bnd_flag,darcy_class)
%% function to compute inflow/outflow of boundary-lying elements

if sum(bnd_flag)==3
    %error('This code does not support an triangular element with three boundary nodes.');
end
qn = zeros(3,1);

o = sum(x)/3;
a = sum(x([1 2],:))/2;
b = sum(x([2 3],:))/2;
c = sum(x([1 3],:))/2;

n1 = [o(2) - a(2), a(1)-o(1)];
n2 = [o(2) - b(2), b(1)-o(1)];
n3 = [o(2) - c(2), c(1)-o(1)];

%% The element contains an edge lying on the inlet
if sum(inlet_flag) == 2 
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
      
    %% At this moment, we ignore elements having single inlet node.
elseif sum(inlet_flag) == 1 
    if sum(1-fFactor < eps)>= 1
        qn(1) = -(n1-n3)*v;
        qn(2) = -(n2-n1)*v;
        qn(3) = -(n3-n2)*v;
    end
else
    error('The triangular element must have at most two nodes on the inlet boundary.');
end

qn = qn*darcy_class.thickness;

end