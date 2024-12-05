function u = sim_permeability(u_bar,C)

u = mvnrnd(u_bar,C,1);

plot(u)
ylim([-2,2])

end