function [p, surfaceCell,faceNormal] = buckyballs()
% Ref. to: https://www.goldennumber.net/bucky-balls/
%% Author: Wan Ji, Wuhan University, Wuhan, China
% Any advice or questions, please contact me via E-mail:
% wanji@whu.edu.cn
% Date of this version 2021/06/28
phi = (sqrt(5)+1)/2;
p = [
0,1,3*phi;
0,1,-3*phi
0,-1,3*phi
0,-1,-3*phi
1,3*phi,0
1,-3*phi,0
-1,3*phi,0
-1,-3*phi,0
3*phi,0,1
3*phi,0,-1
-3*phi,0,1
-3*phi,0,-1
2,(1+2*phi),phi
2,(1+2*phi),-phi
2,-(1+2*phi),phi
2,-(1+2*phi),-phi
-2,(1+2*phi),phi
-2,(1+2*phi),-phi
-2,-(1+2*phi),phi
-2,-(1+2*phi),-phi
(1+2*phi),phi,2
(1+2*phi),phi,-2
(1+2*phi),-phi,2
(1+2*phi),-phi,-2
-(1+2*phi),phi,2
-(1+2*phi),phi,-2
-(1+2*phi),-phi,2
-(1+2*phi),-phi,-2
phi,2,(1+2*phi)
phi,2,-(1+2*phi)
phi,-2,(1+2*phi)
phi,-2,-(1+2*phi)
-phi,2,(1+2*phi)
-phi,2,-(1+2*phi)
-phi,-2,(1+2*phi)
-phi,-2,-(1+2*phi)
1,(2+phi),2*phi
1,(2+phi),-2*phi
1,-(2+phi),2*phi
1,-(2+phi),-2*phi
-1,(2+phi),2*phi
-1,(2+phi),-2*phi
-1,-(2+phi),2*phi
-1,-(2+phi),-2*phi
(2+phi),2*phi,1
(2+phi),2*phi,-1
(2+phi),-2*phi,1
(2+phi),-2*phi,-1
-(2+phi),2*phi,1
-(2+phi),2*phi,-1
-(2+phi),-2*phi,1
-(2+phi),-2*phi,-1
2*phi,1,(2+phi)
2*phi,1,-(2+phi)
2*phi,-1,(2+phi)
2*phi,-1,-(2+phi)
-2*phi,1,(2+phi)
-2*phi,1,-(2+phi)
-2*phi,-1,(2+phi)
-2*phi,-1,-(2+phi)
];
col = jet(200);
col = col(randperm(200),:);
[surfaceCell,faceNormal] = surfaceReduction(p, 1:1:size(p,1),1e-12);

figure(1)
clf
for i = 1:1:numel(surfaceCell)
    k = surfaceCell{i,1};
    patch('vertices',p(k,:),'faces',1:1:numel(k),'FaceColor',col(i,:),'facealpha',0.5)
    hold on
% if(numel(k)==6)
%      patch('vertices',p(k,:),'faces',1:1:numel(k),'FaceColor','k','facealpha',1,'edgecolor','k')
%      hold on
% else
%          patch('vertices',p(k,:),'faces',1:1:numel(k),'FaceColor','w','facealpha',1,'edgecolor','k')
%      hold on
% end
end
view(-67,42);
axis equal
axis off