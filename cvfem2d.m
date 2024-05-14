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

function opt = cvfem2d_init(opt)

opt.cvfem.activeNode = zeros(opt.mesh.nnode,1);
opt.cvfem.activeNode(opt.bndry.inlet_flag==1) = 1;

opt.cvfem.activeElement = zeros(opt.mesh.nelem,1);
inlet_idx = find(opt.bndry.inlet_flag);
% At the begining, there are no active elements because no flow moves into
% the domain through the inlet yet. But in the flux calculation, elements
% involving inlet nodes needs to be highlighted somehow. 
% We assigned 0.5 to these elements to distinguish them from finite elements.
candidate = zeros(opt.mesh.nelem,1);
for i = 1 : nnz(opt.bndry.inlet_flag)
    ival = inlet_idx(i);
    for j = 2 : opt.mesh.has_node_i(ival,1)+1
        candidate(opt.mesh.has_node_i(ival,j))= 1;
    end
end
candidate_idx = find(candidate);
for i = 1 : length(candidate_idx)
    if sum(opt.bndry.inlet_flag(opt.mesh.elem(candidate_idx(i),:)))==2
        opt.cvfem.activeElement(candidate_idx(i)) = 0.5;
    end
end

opt.mesh.elem_including_inlet_edge = sparse(opt.cvfem.activeElement==0.5);

%% filling time and filling factor
opt.cvfem.fTime = 0;                         % The filling time
opt.cvfem.fFactor = zeros(opt.mesh.nnode,1); % the filling factor

%% allocate a sparse matrix 
opt.cvfem.A = spalloc(opt.mesh.nnode,opt.mesh.nnode,opt.mesh.nnode*10);

end