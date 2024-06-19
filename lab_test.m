%% Problem set up
my_mesh = DelaunayMesh(p,e,t);
my_physics = Physics(0.1, 0.35, 0.2, @permeability);
my_pressure = Pressure(my_mesh,my_physics,@is_inlet,@is_vent,@p_D);

%% compile CVFEM method
my_RTMFlow = RTMFlowDryspot(my_mesh,my_physics,my_pressure);

%% Execute solver
my_RTMFlow.run()

%% Argument set up
function K = permeability(x)

    K = 1e-14 * eye(2);

    %% is point in channels
    is_left_channel = abs(x(1)-0.1045)<= 0.002;
    is_right_channel = abs(x(1)-0.3135)<= 0.002;
    is_lower_channel = abs(x(2)-0.1045)<= 0.002;
    is_upper_channel = abs(x(2)-0.3135)<= 0.002;

    if  is_left_channel || is_right_channel || is_upper_channel || is_lower_channel

        K = 1e-10*eye(2);
    end
end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure + 1e5;
p(pressure_class.is_inlet) = 1.5e5;

end

function bool = is_inlet(node)

%% inlets on central circles (eps neeeded for rounding error)
is_north_inlet = norm(node-[0.209, 0.3135]) < 0.0025+eps;
is_south_inlet = norm(node-[0.209, 0.1045]) < 0.0025+eps;
is_east_inlet = norm(node-[0.3135, 0.209]) < 0.0025+eps;
is_west_inlet = norm(node-[0.1045, 0.209]) < 0.0025+eps;

bool = is_north_inlet || is_south_inlet || is_east_inlet || is_west_inlet;

end

function bool = is_vent(node)

bool = (node(1) == 0 || node(2) == 0 || node(1) == 0.418 || node(2) == 0.418);

end