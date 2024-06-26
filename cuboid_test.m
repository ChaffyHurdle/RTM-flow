model = createpde;

%{
Mesh geometries

multicuboid, multicylinder, or multisphere 

mesh = generateMesh(model);
[p,e,t] = meshToPet(mesh); %convert form of mesh

figure; pdemesh(mesh);
figure; pdemesh(p,e,t);

translate(g,[1 2]);

g3 = addCell(g1,g2)
%}

%{
model = createpde;
gm = multisphere(1);
model.Geometry = gm;

gm2 = multicuboid(1,2,1);
translate(gm2,[0.5 1 0]);

model.Geometry = add(gm,gm2);



mesh = generateMesh(model,"Hmin",0.5);
[p,e,t] = meshToPet(mesh);

figure; pdemesh(mesh);

%}

cuboid_mesh

model = createpde;

nodes = msh.POS'; 
groupsID = msh.TETS(:,5);
elements = msh.TETS(:,1:4)';
geom2 = model.geometryFromMesh(nodes,elements); 

%% Plotting:
figure(1); pdeplot3D(model,'FaceAlpha',0.5)                     %Mesh Plot
figure(2); pdegplot(model,'CellLabels','on','FaceAlpha',0.5)    %Geometry
figure(3); pdegplot(model,'FaceLabels','on','FaceAlpha',0.5)    %Geometry

