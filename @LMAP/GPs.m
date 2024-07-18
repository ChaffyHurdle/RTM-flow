% Define the number of points and their coordinates
x = my_mesh.centroids; % Random points in 2D for example
N = length(x); % Number of points

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
        K(i, j) = covariance_function(x(i, :), x(j, :));
    end
end

% Add a small noise term to the diagonal for numerical stability
K = K + 1e-6 * eye(N);

% Sample from the multivariate normal distribution
L = chol(K, 'lower'); % Cholesky decomposition
z = randn(N, 1); % Standard normal random variables
f_sample = m + L * z; % Sample from the GP

% Define the grid
gridSize = 100; % Grid size
xq = linspace(min(x(:,1)), max(x(:,1)), gridSize);
yq = linspace(min(x(:,2)), max(x(:,2)), gridSize);
[Xq, Yq] = meshgrid(xq, yq);

% Interpolate the GP sample to the regular grid
F = scatteredInterpolant(x(:,1), x(:,2), f_sample, 'linear', 'none');
Vq = F(Xq, Yq);

% Plot the heatmap
figure;
imagesc(xq, yq, exp(Vq));
set(gca, 'YDir', 'normal');
colorbar;
title('Heatmap of GP Sample');
xlabel('x');
ylabel('y');

