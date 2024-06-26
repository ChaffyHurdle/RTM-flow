clc;clear
%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
%% Example
Nodes = [
0,0,0;
2,0,0;
1,1.5,0;
1,1,1.5;
1,0,1;
];
Faces = {
[1,2,3];[2,3,4];[3,1,4];[1,2,4];
[1,2,5];[2,4,5];[4,1,5];
};
% patch two tetrahedrons(Of course polyhedron like voronoi can be patched like this)
Voronois = {[1,2,3,4];[4,5,6,7]};
v3dCell = voronoi3d(Nodes, Faces, Voronois);
[c,d] = voronoiPatch(v3dCell,'bycell');
[a,b] = v3dCell.voronoiPatch('nodal',Nodes(:,3));