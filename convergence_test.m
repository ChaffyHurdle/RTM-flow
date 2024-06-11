%% Problem set up
my_mesh = DelaunayMesh(p,e,t);
my_darcy = Physics(0.1, 0.35, 1, @permeability);
my_pressure = Pressure(my_mesh,my_darcy,@is_inlet,@is_vent,@p_D);

%% compile RTM-flow method
my_RTMflow = RTMFlowConvergence(my_mesh,my_darcy,my_pressure);
my_RTMflow.visualise_class.is_plotting_volume = true;

%% Execute solver
[mid_RTMflow, final_RTMflow] = my_RTMflow.run(0.5);

filling_errors(1) = compute_filling_error(final_RTMflow);
[pressure_errors(1), boundary_errors(1)] = compute_flow_error(mid_RTMflow); 

%% Analytical solutions
function error = compute_filling_error(RTMflow_class)

%% extract solution properties
time = RTMflow_class.time;
physics_class = RTMflow_class.physics_class;
pressure_class = RTMflow_class.pressure_class;

K = physics_class.permeability([0 0]);
K = K(1,1);
mu = physics_class.viscosity;
phi = physics_class.porosity;
p_0 = max(pressure_class.pressure);

%% compute filling time error
filling_time = (phi*mu)/(2*K*p_0);
error = abs(filling_time - time)/filling_time;
disp(['filling time error = ' num2str(error)])

end

function [p_error, x_error] = compute_flow_error(RTMflow_class)

time = RTMflow_class.time;
physics_class = RTMflow_class.physics_class;
pressure_class = RTMflow_class.pressure_class;

K = physics_class.permeability([0 0]);
K = K(1,1);
mu = physics_class.viscosity;
phi = physics_class.porosity;
p_0 = max(pressure_class.pressure);

%% compute boundary error
x_front = sqrt((2*K*p_0*time)/(phi*mu));
is_moving_boundary = pressure_class.is_Dirichlet...
                   & pressure_class.is_node_active...
                   & ~pressure_class.is_inlet;
boundary_positions = pressure_class.mesh_class.nodes(is_moving_boundary,:);
x_front_h = boundary_positions(:,1);
x_front_diff = abs(x_front_h - x_front);
x_error = mean(x_front_diff);
disp(['boundary error = ' num2str(x_error)])

%% compute pressure error
x_vec = RTMflow_class.Delaunay_mesh_class.nodes(:,1);
p_true = 0*x_vec;
p_true(x_vec <= x_front) = (p_0 * (1 - x_vec(x_vec <= x_front)./x_front));
p_error = mean(abs(p_true - pressure_class.pressure))/p_0;
disp(['pressure error = ' num2str(p_error)])

end

%% Argument set up
function K = permeability(x)
    
    K = 10e-10 * eye(2);

end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure;
p(pressure_class.is_inlet) = 1.5e5;

end

function bool = is_inlet(node)

bool = (node(1) == 0);

end

function bool = is_vent(node)

bool = (node(1) == 1);

end