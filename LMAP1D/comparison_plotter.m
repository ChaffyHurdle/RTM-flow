function pl = comparison_plotter(u_map,C_map,U_EKI,U_RML,U_MCMC,params,experiment)

hex = "#d3d3d3";

figure(5)

% LMAP plot
subplot(2,2,1)
upper50 = u_map + 0.674*sqrt(diag(C_map))';
lower50 = u_map - 0.674*sqrt(diag(C_map))';
upper95 = u_map + 1.96*sqrt(diag(C_map))';
lower95 = u_map - 1.96*sqrt(diag(C_map))';
fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none')
hold on
plot(params.x_locations,u_map,"k--")
plot(params.x_locations,upper95,"k")
plot(params.x_locations,lower95,"k")
plot(experiment.params_fwd.x_locations,experiment.u_true,"r")
xline(experiment.ups_true(end),"r--")
hold off
xlim([0,1])
ylim([-2,2])
xlabel('$$x$$','interpreter','latex')
ylabel('$$u(x)$$','interpreter','latex')
title('LMAP','interpreter','latex')

% EKI plot
subplot(2,2,2)
U_mean=mean(U_EKI,2);
upper50 = U_mean' + 0.674*sqrt(var(U_EKI,0,2))';
lower50 = U_mean' - 0.674*sqrt(var(U_EKI,0,2))';
upper95 = U_mean' + 1.96*sqrt(var(U_EKI,0,2))';
lower95 = U_mean' - 1.96*sqrt(var(U_EKI,0,2))';
fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none');
hold on
plot(params.x_locations,U_mean,'k--')
plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
plot(params.x_locations,upper95,'k')
plot(params.x_locations,lower95,'k')
xline(experiment.ups_true(end),"r--")
hold off
%legend('$$U_{EKI}\pm 0.674\sigma_{EKI}$$','$$U_{EKI}$$','truth','location','north','fontsize',20,'interpreter','latex')
xlim([0,1])
ylim([-2,2])
xlabel('$$x$$','interpreter','latex')
ylabel('$$u(x)$$','interpreter','latex')
title('EKI','interpreter','latex')


% RML plot
subplot(2,2,3)
U_mean=mean(U_RML,2);
upper50 = U_mean' + 0.674*sqrt(var(U_RML,0,2))';
lower50 = U_mean' - 0.674*sqrt(var(U_RML,0,2))';
upper95 = U_mean' + 1.96*sqrt(var(U_RML,0,2))';
lower95 = U_mean' - 1.96*sqrt(var(U_RML,0,2))';
fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none');
hold on
plot(params.x_locations,U_mean,'k--')
plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
plot(params.x_locations,upper95,'k')
plot(params.x_locations,lower95,'k')
xline(experiment.ups_true(end),"r--")
hold off
xlim([0,1])
ylim([-2,2])
xlabel('$$x$$','interpreter','latex')
ylabel('$$u(x)$$','interpreter','latex')
title('RML','interpreter','latex')


% MCMC plot
subplot(2,2,4)
U_mean=mean(U_MCMC,2);
upper50 = U_mean' + 0.674*sqrt(var(U_MCMC,0,2))';
lower50 = U_mean' - 0.674*sqrt(var(U_MCMC,0,2))';
upper95 = U_mean' + 1.96*sqrt(var(U_MCMC,0,2))';
lower95 = U_mean' - 1.96*sqrt(var(U_MCMC,0,2))';
fill([params.x_locations,fliplr(params.x_locations)], [lower50,fliplr(upper50)], hex2rgb(hex), 'EdgeColor','none');
hold on
plot(params.x_locations,U_mean,'k--')
plot(experiment.params_fwd.x_locations,experiment.u_true,'r')
plot(params.x_locations,upper95,'k')
plot(params.x_locations,lower95,'k')
xline(experiment.ups_true(end),"r--")
xlabel('$$x$$','interpreter','latex')
ylabel('$$u(x)$$','interpreter','latex')
hold off
xlim([0,1])
ylim([-2,2])
title('MCMC','interpreter','latex')

figure(8)
% Means comparison
subplot(1,2,1)
plot(params.x_locations, u_map,"red")
hold on
plot(params.x_locations, mean(U_EKI,2),"blue")
plot(params.x_locations, mean(U_RML,2),"green")
plot(params.x_locations, mean(U_MCMC,2),"k--")
hold off
xlim([0,1])
ylim([-2,2])
legend("LMAP","EKI","RML","pcn-MCMC")
title("Means")


% Std comparison
subplot(1,2,2)
plot(params.x_locations, sqrt(diag(C_map)),"red")
hold on
plot(params.x_locations, sqrt(var(U_EKI,0,2)),"blue")
plot(params.x_locations, sqrt(var(U_RML,0,2)),"green")
plot(params.x_locations, sqrt(var(U_MCMC,0,2)),"k--")
hold off
xlim([0,1])
ylim([0,1])
legend("LMAP","EKI","RML","pcn-MCMC")
title("Stds")

end