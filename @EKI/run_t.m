function obj = run_t(obj,t)

%% Runs EKI up to time t \in [1,2,3,4,5]

% Various shorthands
u = obj.inverse_class.u0;
N_En = obj.J;
U = mvnrnd(u,obj.inverse_class.C0_inv,N_En)';
physics_class = obj.physics_class;
mesh_class = obj.mesh_class;

Sigma = reshape(obj.inverse_class.Sigma(:,1:t),[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
d = reshape(obj.inverse_class.data(:,1:t),[],1);

M = length(d);
U_mean = mean(U,2);

t_vec(1) = 0;
Cond = 1;
MAX = obj.max_iterations;
iter = 0;
tic;

% Start EKI
while (Cond==1)&&(iter<MAX)
    iter=iter+1;
    Z = zeros(M,N_En);
    parfor en = 1:N_En
        disp(en)
        K_true = exp(U(:,en));
        physics_class_en = physics_class;
        physics_class_en.permeability = K_true;
        pressure_class_en = Pressure(mesh_class,physics_class_en);
        RTMflow_class_en = RTMFlow(mesh_class,physics_class_en,pressure_class_en,0);
        RTMflow_class_en = RTMflow_class_en.run(physics_class.observation_times(t));
        p = reshape(RTMflow_class_en.pressure_data(:,1:t),[],1);
    
        Z(:,en)=Sigma_minus_half*(d-p);
    end
    Z_m=mean(Z,2);
    Misfit(iter)=norm(Z_m(:,1))^2/M;
    alpha=mean(vecnorm(Z).^2)/M;
    if (t_vec(iter)+1/alpha>1)
        alpha=1/(1-t_vec(iter));
        Cond=0;
    end
    
    Delta_Z=Z-Z_m;
    C_tilde=1/(N_En-1)* (Delta_Z)*(Delta_Z)'+alpha*eye(M);
    Delta_RN=U-U_mean;
    C_RN=1/(N_En-1)*Delta_RN*Delta_Z';
    E=sqrt(alpha)*randn(M,N_En);
    E=E-mean(E,2);
    W=Z+E;
    U=U-C_RN*(C_tilde\W);
    U_mean=mean(U,2);
    t_vec(iter+1)=t_vec(iter)+1/alpha;
    
    % Plot iterate
    figure(1)
    subplot(1,3,1)
    pdeplot(obj.inverse_class.fwd_mesh.nodes',obj.inverse_class.fwd_mesh.elements',XYData=obj.inverse_class.u_true, ...
        XYStyle='interp',ColorMap="turbo",Mesh="off")
    clim([-1.5,1.5])
    subplot(1,3,2)
    pdeplot(mesh_class.nodes',mesh_class.elements',XYData=U_mean, ...
        XYStyle='interp',ColorMap="turbo",Mesh="off")
    clim([-1.5,1.5])
    subplot(1,3,3)
    pdeplot(mesh_class.nodes',mesh_class.elements',XYData=var(U,0,2), ...
        XYStyle='interp',ColorMap="turbo",Mesh="off")
    clim([0,0.25])
    drawnow
    
end

% Save data
obj.F_evals = iter*N_En;
timer = toc
obj.timer = timer;
obj.ueki = U_mean;
obj.Ceki = var(U,0,2);
end
