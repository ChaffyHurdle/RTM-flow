%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
clc;clear
%% sub Exam01
p = icosohedron();
%% sub Exam02
% p = rand(50,3);
%% get the delaunayTriangulation geometry
a = delaunayTriangulation(p);
elements = a.ConnectivityList';
nodes = a.Points';
model = createpde();
geometryFromMesh(model,nodes,elements);
pdegplot(model,'FaceLabels','on','FaceAlpha',0.5)
generateMesh(model,'hmax',0.5);
% pdeplot3D(model)

obj = voronoi3d; % create an empty voronoi3d object
obj = createVoronoi3dFromMesh(obj, model.Mesh.Nodes', model.Mesh.Elements'); % create voronoi by geometry
[c,d] = voronoiPatch(obj,'bycell');
hold on
scatter3(obj.Nodes(:,1),obj.Nodes(:,2),obj.Nodes(:,3),10,'k','filled')
scatter3(obj.Centroids(:,1),obj.Centroids(:,2),obj.Centroids(:,3),20,'r','filled')
axis off

[c,d] = voronoiPatch(obj,'nodal',obj.Nodes(:,2));