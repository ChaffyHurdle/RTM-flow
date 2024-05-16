function g = g_D(points,bndry)

g = zeros(size(points,1),1);
g(points(:,2)==0) = 1.5e5;

end