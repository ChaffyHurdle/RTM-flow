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
for i = 1:4

    %% set up
    ram = whos;
    initial_RAM_used = sum([ram.bytes])/1000000;
    tic

    mesh = DelaunayMesh(mesh_files(i));
    pressure = Pressure(mesh,darcy_rules,@is_inlet,@is_vent,@p_D);

    %% compile RTM-flow method
    RTM = RTMFlow(mesh,darcy_rules,pressure);

    %% run simulation
    RTM = RTM.run();

    %% store results
    run_time(i) = toc;
    ram = whos;
    final_RAM_used = sum([ram.bytes])/1000000;
    memory_used(i) = final_RAM_used - initial_RAM_used;

    clear RTM
    clear mesh
    clear pressure


end

T_1 = table(num_dofs,run_time,memory_used)
writetable(T_1,'unit_sqr_performance.csv')

%% Run test
for i = 1:4

    %% set up
    ram = whos;
    initial_RAM_used = sum([ram.bytes])/1000000;
    tic

    mesh = DelaunayMesh(mesh_files(i));
    pressure = Pressure(mesh,darcy_rules,@is_inlet,@is_vent,@p_D);

    %% compile RTM-flow method
    RTM = RTMFlow(mesh,darcy_rules,pressure);
    RTM.visualise_class.is_plotting_volume = true;

    %% run simulation
    RTM = RTM.run();

    %% store results
    run_time(i) = toc;
    ram = whos;
    final_RAM_used = sum([ram.bytes])/1000000;
    memory_used(i) = final_RAM_used - initial_RAM_used;

    clear RTM
    clear mesh
    clear pressure


end

T_2 = table(num_dofs,run_time,memory_used)
writetable(T_2,'unit_sqr_plot_performance.csv')



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