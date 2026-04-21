function C = compute_matern_covariance(obj,centroids)

% Constructs Matern covariance matrix

xxyy = centroids;
N = length(xxyy);

% Define the covariance function
sigma_f = obj.matern_var;    % Variance
l = obj.matern_length_scale; % Length scale
nu = obj.matern_nu;          % Smoothness
covariance_function = @(x1, x2) sigma_f * (2^(1-nu)/gamma(nu)) * ((sqrt(2*nu)*norm(x1 - x2)/l)^nu) * besselk(nu,sqrt(2*nu)*norm(x1 - x2)/l);

% Compute the covariance matrix C
C = zeros(N, N);
for i = 1:N
    for j = (i+1):N
        C(i,j) = covariance_function(xxyy(i,:), xxyy(j,:));
        C(j,i) = C(i,j);
    end
    C(i,i) = sigma_f;
end


