clc;clear
%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
%% Examples
% myGeo3d = multisphere(1); % create a sphere geometry with radius 1
 myGeo3d = multicylinder(1,3); % create a cylinder geometry with radius 1 and height 3
% myGeo3d = multicuboid(1,1,1); % create a cuboid geometry 1*1*1
obj = voronoi3d; % create an empty voronoi3d object
obj = createVoronoi3dFromGeometry(obj, myGeo3d, 0.3); % create voronoi by geometry
voronoiPatch(obj,'bycell');
hold on
scatter3(obj.Nodes(:,1),obj.Nodes(:,2),obj.Nodes(:,3),10,'k','filled')
scatter3(obj.Centroids(:,1),obj.Centroids(:,2),obj.Centroids(:,3),20,'r','filled')
axis off

voronoiPatch(obj,'nodal',obj.Nodes(:,2));

% q = arrayfun(@(i)numel(obj.Voronois{i}),(1:1:numel(obj.Voronois))');
% cellNum = find(q==max(q));
% Faces = obj.Voronois{cellNum};
% FNormal =  obj.FaceNormals(Faces,:);