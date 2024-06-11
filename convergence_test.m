%% Problem set up
my_mesh = DelaunayMesh(p,e,t);
my_darcy = Physics(0.1, 0.35, 1, @permeability);
my_pressure = Pressure(my_mesh,my_darcy,@is_inlet,@is_vent,@p_D);

%% compile RTM-flow method
my_RTMflow = RTMFlow(my_mesh,my_darcy,my_pressure);

%% Execute solver
my_RTMflow = my_RTMflow.run();

compute_error(my_RTMflow)

%% Analytical solutions
function compute_error(RTMflow_class)

time = RTMflow_class.time;
physics_class = RTMflow_class.physics_class;
pressure_class = RTMflow_class.pressure_class;

K = physics_class.permeability([0 0]);
K = K(1,1);
mu = physics_class.viscosity;
phi = physics_class.porosity;
p_0 = max(pressure_class.pressure);

filling_time = (phi*mu)/(2*K*p_0);

filling_time_error = abs(filling_time - time)/filling_time;

x_front = sqrt((2*K*p_0*filling_time)/(phi*mu));

is_moving_boundary = pressure_class.is_Dirichlet...
                   & pressure_class.is_node_active...
                   & ~pressure_class.is_inlet;
boundary_positions = pressure_class.mesh_class.nodes(is_moving_boundary,:);
x_front_h = boundary_positions(:,1);

x_front_diff = abs(x_front_h - x_front);
x_front_error = mean(x_front_diff);

disp(['boundary error = ' num2str(x_front_error)])
disp(['filling time error = ' num2str(filling_time_error)])

end

%% Argument set up
function K = permeability(x)
    
    K = 10e-10 * eye(2);

end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure;
p(pressure_class.is_inlet) = 1;

end

function bool = is_inlet(node)

bool = (node(1) == 0);

end

function bool = is_vent(node)

bool = (node(1) == 1);

end