function C0 = compute_matern_covariance(obj,centroids)

% Define the number of points on highly refined structured 2D grid
% x = linspace(0,1,obj.high_ref);
% y = x;
% [xx,yy] = meshgrid(x,y);
% xx = reshape(xx,[],1);
% yy = reshape(yy,[],1);
% xxyy = [xx yy];
xxyy = centroids;
N = length(xxyy);

% Define the covariance function (Squared Exponential/RBF Kernel)
sigma_f = obj.matern_var; % Signal variance
l = obj.matern_length_scale; % Length scale
nu = obj.matern_nu;

covariance_function = @(x1, x2) sigma_f * (2^(1-nu)/gamma(nu)) * ((sqrt(2*nu)*norm(x1 - x2)/l)^nu) * besselk(nu,sqrt(2*nu)*norm(x1 - x2)/l);

% Compute the covariance matrix K
C = zeros(N, N);
for i = 1:N
    for j = (i+1):N
        C(i,j) = covariance_function(xxyy(i,:), xxyy(j,:));
        C(j,i) = C(i,j);
    end
    C(i,i) = sigma_f;
end
C0 = C;


