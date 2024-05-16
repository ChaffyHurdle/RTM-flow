 close all
% Import mesh files
my_mesh = Mesh(p,e,t);

% Darcy rules set
my_darcy = Darcy(0.1, 0.35, 0.2, @permeability);

my_pressure = Pressure(my_mesh,'inlet_location','vent_location',@p_D);
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

function p = p_D(points,bndry,time)

p = zeros(size(points,1),1);
p(points(:,2)==0) = 1.5e5;

end
    
    





    