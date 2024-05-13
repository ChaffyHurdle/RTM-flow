close all
% Import mesh files
my_mesh = Mesh(p,e,t);

% Darcy rules set
my_darcy = Darcy(0.1, 0.35, 0.2, @permeability);


% Boundary conditions set
bndry.gd_filename   = 'g_D';
bndry.gn_filename   = '';
bndry.pvent         = 0;
bndry.pinlet        = 1.5e5;
bndry.inlet_location_fname = 'inlet_location';
bndry.vent_location_fname = 'vent_location';

% First set up routine
opt = cvfem_setup(mesh,bndry,darcy,@permeability);

% Execute code
opt = cvfem2d(opt);

%% Argument set up
function K = permeability(x)
    
    K = [1e-10 0; 0 1e-10];
end

function pressure_class = boundary_pressure(mesh_class,pressure_class)

    nodes = mesh_class.nodes; time = pressure_class.time;
    inlet_nodes = pressure_class.inlet_nodes;

    pressure_class.pressure(:) = 0;

    for i = inlet_nodes
        pressure_class.pressure(i) = 1.5e5;
    end

end
    
    





    