function [surfaceCell,faceNormal] = surfaceReduction(nodeIn, faceCells, tol)
%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
% reduce the surface of 3D triangle to 3D polygon
% sometimes like a co_planar use
% input variables: 
% nodeIn: the node for all the voronoi cell
% faceCells:  a m-by-3 node-ID matrix that forms the voronoi cell
% tol: tolerance for coplanar judge
% outPut variables:
% surfaceCell: the Final surface that diminishes the coplanar triangle
% faceNormal: normal of each 3d polygon
n = faceCells;
p = nodeIn(n,:);
[k,~] = convhulln(p);
mp = n(k);
normalDir  = zeros(size(k,1),3);
[~, ~, Cx, Cy, Cz]=Polyhedron_VAC(p,k);
centro = [Cx,Cy,Cz];
for j = 1:1:size(k,1)
    tri = k(j,:);
    x = p(tri,1);
    y = p(tri,2);
    z = p(tri,3);
    p1 = [x(1),y(1),z(1)];
    p2 = [x(2),y(2),z(2)];
    p3 = [x(3),y(3),z(3)];
    normalDir(j,:) = cross(p2-p1,p3-p1);
    normalDir(j,:) = normalDir(j,:)/norm(normalDir(j,:));
    if(([x(1),y(1),z(1)]-centro)*normalDir(j,:)'<0)
        normalDir(j,:) = -normalDir(j,:);
    end
end
[a,~,ic] = uniquetol(normalDir,tol,'ByRows',true); % find the triangles at the same plane
faceNormal = a;
faces = cell(size(a,1),1);
for j = 1:1:numel(ic)
    faces{ic(j),1} = [faces{ic(j),1}; mp(j,:)'];
end
for j = 1:1:numel(faces)
    faces{j,1} = unique(faces{j,1});
end
t = [2,3,1,2];
for j = 1:1:numel(faces)
    p = faces{j,1};
    enode = nodeIn(p,:);
    mDim = find(max(abs(a(j,:)))==abs(a(j,:)));
    t0 = t(mDim(1):mDim(1)+1);
    pp = enode(:,t0);
    kk= convhull(pp);
    faces{j,1} = p(kk);
end
surfaceCell = faces;
end