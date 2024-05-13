function cvFEM_plot(opt)
%This function is designed to produce an animated/still image within MATLAB
%of the cvFEM for each time-step of the simulation

% temporary plotting variables
p = opt.mesh.node';
t = opt.mesh.elem';
pressure = opt.cvfem.u;

% solution plotting
figure(1)
clf
axis equal
axis auto
pdeplot(p,t,XYData=pressure,ZData=pressure,ColorMap="jet")
title(["time = " num2str(opt.cvfem.fTime)])
%pause(1e-10)

% Add axis and title labels
%
%
%

% Add volume density tracking on the 2D mesh
fluid_elements = zeros(size(t,2),1);
fluid_elements(opt.cvfem.activeElement ~= 0) = 1;
figure(2)
clf;
%trisurf(t,p(1,:),p(2,:),fluid_elements)
pdeplot(p,t,XYData=fluid_elements,ColorMap="jet",Mesh="on")
axis equal
title(["time = " num2str(opt.cvfem.fTime)])
%pause(1e-10)
