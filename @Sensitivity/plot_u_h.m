function obj = plot_u_h(obj)

cmin = min([min(obj.u),min(obj.h),min(obj.u + obj.h)]);
cmax = max([max(obj.u),max(obj.h),max(obj.u + obj.h)]);

figure
subplot(1,3,1)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',...
    XYData=obj.u,XYStyle='interp',ColorMap="jet",Mesh="off")
caxis([cmin,cmax])

subplot(1,3,2)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',...
    XYData=obj.h,XYStyle='interp',ColorMap="jet",Mesh="off")

subplot(1,3,3)
pdeplot(obj.mesh_class.nodes',obj.mesh_class.elements',...
    XYData=obj.u+obj.h,XYStyle='interp',ColorMap="jet",Mesh="off")
caxis([cmin,cmax])


end