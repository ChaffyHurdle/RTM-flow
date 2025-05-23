function [U,timer,iter] = run_EKI(u0,C,params,experiment,N_En,plotting)

C_minus_half = inv(sqrtm(C));
u_list = u0;
U = mvnrnd(u0,C,N_En)';

Sigma = reshape(experiment.Sigma,[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
Sigma = diag(Sigma);
d = reshape(experiment.d,[],1);

M=length(d);


delete(gcp('nocreate'))
i=10;
fprintf('Number of slots available: %d\n', i);
parpool('Threads', i);

U_mean=mean(U,2);

t(1)=0;
Cond=1;
MAX=120;
iter=0;
tic;
while (Cond==1)&&(iter<MAX)
    iter=iter+1
    Z=zeros(M,N_En);
    parfor en=1:N_En
        [p,~] = forward_map(U(:,en)',params);
        Z(:,en)=Sigma_minus_half*(d-reshape(p,[],1));
    end
    Z_m=mean(Z,2);
    Misfit(iter)=norm(Z_m(:,1))^2/M;
    alpha=mean(vecnorm(Z).^2)/M;
    if (t(iter)+1/alpha>1)
        alpha=1/(1-t(iter));
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
    t(iter+1)=t(iter)+1/alpha;
    if plotting
        figure(2)
        X=[params.x_locations,fliplr(params.x_locations)];
        Y=[(U_mean-0.674*sqrt(var(U,0,2)))',fliplr((U_mean+0.674*sqrt(var(U,0,2)))')];
        fill(X,Y, 'g', 'EdgeColor','none', 'FaceAlpha',0.25);
        hold on
        plot(params.x_locations,U_mean,'k--')
        plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
        plot(params.x_locations,U_mean+1.96*sqrt(var(U,0,2)),'k')
        plot(params.x_locations,U_mean-1.96*sqrt(var(U,0,2)),'k')
        xline(experiment.ups_true(end),"r--")
        %legend('$$U_{EKI}\pm 0.674\sigma_{EKI}$$','$$U_{EKI}$$','truth','location','north','fontsize',20,'interpreter','latex')
        xlim([0,1])
        ylim([-2,2])
        hold off
        drawnow
    end
    
end
timer = toc
iter = iter * N_En;

poolobj = gcp('nocreate');
delete(poolobj);
disp ' ... all done.'
end
