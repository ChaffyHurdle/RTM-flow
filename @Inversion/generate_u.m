function obj = generate_u(obj)
    
% Sample from the multivariate normal distribution
x = linspace(0,1,obj.high_ref);
y = x;
Nx = length(x);

L = obj.cholesky_L; % Cholesky decomposition
NL = length(L);
z = randn(NL, 1); % Standard normal random variables
f_sample = L * z; % Sample from the GP (zero-mean)
prior_highdim = reshape(f_sample,Nx,Nx);

% Interpolate the GP sample to the regular grid
prior_mesh_centers = interp2(x,y,prior_highdim,obj.DelaunayMesh.centroids(:,1),obj.DelaunayMesh.centroids(:,2),'nearest');

obj.u_highdim = prior_highdim;
obj.u_meshcenters = prior_mesh_centers;

% Plot the heatmap
figure;
subplot(1,2,1)
imagesc(prior_highdim);
colormap jet
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
subplot(1,2,2)
scatter(obj.DelaunayMesh.centroids(:,1),obj.DelaunayMesh.centroids(:,2),[],prior_mesh_centers,'filled');
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');

end