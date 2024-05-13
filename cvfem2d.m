function opt = cvfem2d(opt)

%% sets users made file to function for dirichlet boundary
gD = str2func(opt.bndry.gd_filename);

% pressure stored as u???
opt.cvfem.u = zeros(opt.mesh.nnode,1);

%% sets Dirichlet conditions on pressure on inlet
opt.cvfem.u(opt.bndry.inlet_flag==1) = gD(opt.bndry.inlet_pos,opt.bndry);

%% sets pressure on rest of domain to vent pressure
opt.cvfem.u(opt.cvfem.u==0) = opt.bndry.pvent;

%% Initialise activeNode, activeElement vectors, fTime, fFactor, and
%% allocate matrices
opt = cvfem2d_init(opt);

%% Find initial flows out of elements
Q   = update_flow_rate_tri(opt);

%% First numerical time step 
dt  = compute_time_increment(Q,opt.cvfem.fFactor,opt.cvfem.V);
opt.cvfem.fTime = opt.cvfem.fTime+dt;

%% Compute new flow volumes and moving boundaries
opt = update_filling_factor(opt,Q,dt);

cvFEM_plot(opt)
 
while ~isFilled(opt.cvfem.fFactor,opt.mesh.nnode,opt.bndry.vent_idx)

    %% Solve Pressure Problem
    [opt.cvfem, opt.bndry] = still_solver(opt.cvfem,opt.mesh,opt.bndry);
    
    %% Solve flow problem
    Q = update_flow_rate_tri(opt);
    
    %% Compute new time step
    dt = compute_time_increment(Q,opt.cvfem.fFactor,opt.cvfem.V);
    
    %% If nothing to be filled end function
    if dt == 0
        disp('ending due to no boundary motion')
        break;
    end

    %% New time level
    opt.cvfem.fTime = opt.cvfem.fTime + dt;

    %% Update flow volumes and moving boundaries
    opt = update_filling_factor(opt,Q,dt);

    %% Drawing
    cvFEM_plot(opt)
end

if isFilled(opt.cvfem.fFactor,opt.mesh.nnode,opt.bndry.vent_idx)
    disp('ending due to filled mold')
    cvFEM_plot(opt)
end

disp("end")

end


%% This function returns true when each volume is sufficiently full
function flag = isFilled(fFactor,N,vent_idx)
flag = all(abs(fFactor(vent_idx)-1) < eps);
end

function dt = compute_time_increment(Q,fFactor,V)
nfFactor = 1 - fFactor;

candidates = find((Q>eps) & (nfFactor > eps));
if isempty(candidates)
    warning('No positive flux found.')
    dt = 0;
else
    dt = min(nfFactor(candidates).*V(candidates)./Q(candidates));
end
end