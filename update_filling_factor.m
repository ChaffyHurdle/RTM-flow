function opt = update_filling_factor(opt,Q,dt)

%% H.W: unpack variables for now
fFactor = opt.cvfem.fFactor;
V = opt.cvfem.V;

%% increase volume fill factor for any inflow elements
pos_q_idx = find(Q>0);

fFactor_new = min(ones(length(pos_q_idx),1),fFactor(pos_q_idx)+ dt*Q(pos_q_idx)./V(pos_q_idx));

fFactor_new(1-fFactor_new<1e+3*eps) = 1;

new_filled_volume = pos_q_idx(fFactor_new>fFactor(pos_q_idx)& fFactor_new == 1);

fFactor(pos_q_idx) = fFactor_new;
%% H.W: repacking final outputs for now
opt.cvfem.fFactor = fFactor;
opt.cvfem.new_filled_volume = new_filled_volume;