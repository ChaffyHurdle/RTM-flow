function obj = solve(obj)

%% Unpacking all classses

%% Find initial flows out of elements
obj = obj.compute_flow_rates();

%% First numerical time step 
obj = obj.compute_time_increment();
obj.time = obj.time + obj.time_step;

%% Compute new flow volumes and moving boundaries
obj = obj.update_filling_percentage();

while ~isFilled(opt.cvfem.fFactor,opt.mesh.nnode,opt.bndry.vent_idx)

    %% Solve Pressure Problem
    [cvfem, bndry] = update_comp_domain(cvfem,mesh,bndry);
    [opt.cvfem, opt.bndry] = still_solver(opt.cvfem,opt.mesh,opt.bndry);
    
    %% Solve flow problem
    obj = obj.compute_flow_rates();
    
    %% Increment to new time
    obj = obj.compute_time_increment();
    obj.time = obj.time + obj.time_step;

    %% Update flow volumes and moving boundaries
    obj = obj.update_filling_percentage();

end

disp("end")

end

function flag = isFilled(fFactor,N,vent_idx)
flag = all(abs(fFactor(vent_idx)-1) < eps);
end