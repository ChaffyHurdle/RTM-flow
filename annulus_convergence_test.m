%% Problem set up
darcy_rules = Physics(0.1, 0.35, 1, @permeability);

%% Error set up

%% Test mesh strings
mesh_files = ["annulus_mesh_132.mat",...
              "annulus_mesh_478.mat",...
              "annulus_mesh_1812.mat",...
              "annulus_mesh_7048.mat"];
num_dofs = [177, 665, 2577, 10145]';

filling_errors = zeros(4,1);
pressure_errors = zeros(4,1);
boundary_errors = zeros(4,1);

%% Run test
for i = 1:2

    mesh = DelaunayMesh(mesh_files(i));
    pressure = Pressure(mesh,darcy_rules,@is_inlet,@is_vent,@p_D);
    %% compile RTM-flow method
    RTM = RTMFlowConvergence(mesh,darcy_rules,pressure);
    %RTM.visualise_class.is_plotting_volume=true;

    %% run simulation
    [mid_RTMflow, final_RTMflow] = RTM.run(100);

    %% compute errors
    filling_errors(i) = compute_filling_error(final_RTMflow);
    [pressure_errors(i), boundary_errors(i)] = compute_flow_error(mid_RTMflow);

end

T = table(num_dofs,filling_errors,pressure_errors,boundary_errors)
writetable(T,'annulus_error.csv')

%% Analytical solutions
function error = compute_filling_error(RTMflow_class)

%% extract solution properties
time = RTMflow_class.time
physics_class = RTMflow_class.physics_class;
pressure_class = RTMflow_class.pressure_class;

K = physics_class.permeability([0 0]);
K = K(1,1);
mu = physics_class.viscosity;
phi = physics_class.porosity;
p_0 = max(pressure_class.pressure);

%% compute filling time error
filling_time = (phi*mu*0.5^2)/(2*K*p_0)
error = abs(filling_time - time)/filling_time;

end

function [p_error, r_error] = compute_flow_error(RTMflow_class)

time = RTMflow_class.time;
physics_class = RTMflow_class.physics_class;
pressure_class = RTMflow_class.pressure_class;

K = 1e-10;
mu = physics_class.viscosity;
phi = physics_class.porosity;
p_0 = max(pressure_class.pressure);

%% compute boundary error
r_front = sqrt((2*K*p_0*time)/(phi*mu)) + 0.5;
is_moving_boundary = pressure_class.is_Dirichlet...
                   & pressure_class.is_node_active...
                   & ~pressure_class.is_inlet;
boundary_positions = pressure_class.mesh_class.nodes(is_moving_boundary,:);
r_front_h = vecnorm(boundary_positions,2,2);
r_front_diff = abs(r_front_h - r_front);
r_error = mean(r_front_diff);

%% compute pressure error
r_vec = vecnorm(RTMflow_class.Delaunay_mesh_class.nodes,2,2);
p_true = 0*r_vec;
p_true(r_vec <= r_front) = (p_0 * (1 - (r_vec(r_vec <= r_front)-0.5)./(r_front-0.5)));
p_error = mean(abs(p_true - pressure_class.pressure))/p_0;

end

%% Argument set up
function K = permeability(x)

    K = 1e-10*eye(2);

end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure;
p(pressure_class.is_inlet) = 1.5e5;

end

function bool = is_inlet(node)

bool = (norm(node) <= 0.5+eps);

end

function bool = is_vent(node)

bool = (norm(node) >= 1-eps);

end