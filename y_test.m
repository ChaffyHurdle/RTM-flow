%% Problem set up
my_mesh = Mesh(p,e,t);
my_darcy = Darcy(0.1, 0.35, 0.2, @permeability);
my_pressure = Pressure(my_mesh,my_darcy,@is_inlet,@is_vent,@p_D);
my_volume = Volume(my_mesh,my_darcy);
my_velocity = Velocity(my_volume,my_darcy);
my_visuals = Visualisation();

%% compile CVFEM method
my_cvfem = CVFEM(my_mesh,my_pressure,my_volume,my_velocity,my_darcy,my_visuals,[]);

%% Execute solver
my_cvfem.run()

%% Argument set up
function K = permeability(x)
    
    if x(1) < 0
        K = 1e-12*eye(2);
    else
        K = 1e-8*eye(2);
    end
end

function p = p_D(pressure_class)

p = 0*pressure_class.pressure;
p(pressure_class.inlet_flag) = 1.5e5;

end

function bool = is_inlet(node)

bool = (node(2) == 1);

end

function bool = is_vent(node)

bool = (node(2) == -0.9);

end