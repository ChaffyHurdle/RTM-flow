%% Problem set up
clear all
darcy_rules = Physics(0.1, 0.35, 1, @permeability);

%% Error set up

%% Test mesh strings
mesh_files = ["unit_sqr_mesh_177.mat",...
              "unit_sqr_mesh_665.mat",...
              "unit_sqr_mesh_2577.mat",...
              "unit_sqr_mesh_10145.mat"];
num_dofs = [177, 665, 2577, 10145]';

run_time = zeros(4,1);
memory_used = zeros(4,1);
run_time_with_plot = zeros(4,1);
memory_used_with_plot = zeros(4,1);

%% Run test
for i = 1:2

    %% set up
    initial_memory = memory;
    initial_RAM_used = initial_memory.MemUsedMATLAB;
    tic

    mesh = DelaunayMesh(mesh_files(i));
    pressure = Pressure(mesh,darcy_rules,@is_inlet,@is_vent,@p_D);

    %% compile RTM-flow method
    RTM = RTMFlow(mesh,darcy_rules,pressure);

    %% run simulation
    RTM = RTM.run();

    %% store results
    run_time(i) = toc;
    final_memory = memory;
    final_RAM_used = final_memory.MemUsedMATLAB;
    memory_used(i) = final_RAM_used - initial_RAM_used;

    clear RTM
    clear mesh
    clear pressure


end

T = table(num_dofs,run_time,memory_used)
writetable(T,'unit_sqr_performance.csv')


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