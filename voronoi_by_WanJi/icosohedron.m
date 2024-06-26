function [p, surfaceCell,faceNormal] = icosohedron(r)
%(±m,0,±n), (0,±n,±m), (±n,±m,0)
%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
m = sqrt(50-10*sqrt(5))/10;
n = sqrt(50+10*sqrt(5))/10;
p = [
   -n,-m,0; 
   -n,m,0;
   n,-m,0;
   n,m,0;
   0, -n,-m;
   0, -n,m;
   0, n,-m;
   0, n,m;
   -m,0,-n;
   -m,0,n;
   m,0,-n;
   m,0,n+eps;
];
if(nargin==1)
    p = p*r;
end
% [surfaceCell,faceNormal] = surfaceReduction(p, 1:1:size(p,1));
col = jet(20);
col = col(randperm(20),:);
[surfaceCell,faceNormal] = surfaceReduction(p, 1:1:size(p,1), 1e-12);

figure(1)
clf
for i = 1:1:numel(surfaceCell)
    k = surfaceCell{i,1};
    patch('vertices',p(k,:),'faces',1:1:numel(k),'FaceColor',col(i,:),'facealpha',0.5)
    hold on
end
view(-67,42);
axis equal
axis off


