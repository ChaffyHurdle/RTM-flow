function obj = generate_u(obj)

L = obj.cholesky_L; % Cholesky decomposition
NL = length(L);
z = randn(NL, 1); % Standard normal random variables
f_sample = L * z; % Sample from the GP (zero-mean)

obj.u_meshcenters = f_sample;

% Plot the heatmap
figure;
pdeplot(obj.DelaunayMesh.nodes',...
        obj.DelaunayMesh.elements', ...
        XYData=f_sample,XYStyle="flat",ColorMap="jet",Mesh="off")
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');

end