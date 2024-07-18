function [prior_highdim, prior_mesh_centers] = generate_K(ref,my_mesh)

% Define the number of points and their coordinates
x = linspace(0,1,ref); % Random points in 2D for example
y = x;
[xx,yy] = meshgrid(x,y);
xx = reshape(xx,[],1);
yy = reshape(yy,[],1);
xxyy = [xx yy];
N = length(xxyy); % Number of points

% Define the mean function (zero mean in this case)
m = zeros(N, 1);

% Define the covariance function (Squared Exponential/RBF Kernel)
sigma_f = 0.5; % Signal variance
l = 0.1; % Length scale

covariance_function = @(x1, x2) sigma_f^2 * exp(-norm(x1 - x2)^2 / (2 * l^2));

% Compute the covariance matrix K
K = zeros(N, N);
for i = 1:N
    for j = 1:N
        K(i, j) = covariance_function(xxyy(i,:), xxyy(j,:));
    end
end

% Add a small noise term to the diagonal for numerical stability
K = K + 1e-6 * eye(N);

% Sample from the multivariate normal distribution
L = chol(K, 'lower'); % Cholesky decomposition
z = randn(N, 1); % Standard normal random variables
f_sample = m + L * z; % Sample from the GP
prior_highdim = exp(reshape(f_sample,ref,ref));

% Interpolate the GP sample to the regular grid
prior_mesh_centers = interp2(x,y,prior_highdim,my_mesh.centroids(:,1),my_mesh.centroids(:,2),'nearest');

% Plot the heatmap
figure(1);
subplot(1,2,1)
imagesc(prior_highdim);
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');
subplot(1,2,2)
scatter(my_mesh.centroids(:,1),my_mesh.centroids(:,2),[],prior_mesh_centers,'filled');
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');


