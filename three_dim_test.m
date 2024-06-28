clear all; close all; clc;

%% 3D testing

%% set up mesh, physics and pressure
my_mesh = DelaunayMesh3D("one_body.m",10);
my_mesh.plot_mesh();
my_physics = Physics(0.1, 0.35, 1, @permeability);
my_pressure = Pressure3D(my_mesh,my_physics,@is_inlet,@is_vent,@p_D);

%% compile solver
my_RTM = RTMFlow3D(my_mesh,my_physics,my_pressure);

%% Execute solver
my_RTM.run()

%% Argument set up
function K = permeability(x)
    K = 1e-10 * eye(3);
end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure + 1e5;
p(pressure_class.is_inlet) = 1.5e5;

end

function bool = is_inlet(node)

bool = (node(2) >= 49.5);

end

function bool = is_vent(node)

bool = (node(2) == 0);

end