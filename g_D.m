function g = g_D(points,bndry)
%Simple Dirichlet condition of 0 pressure on whole boundary besides the top
%when y == 1, in which case it is fixed to inlet_pressure.

g = zeros(size(points,1),1);
g(points(:,2)==0) = bndry.pinlet;

end