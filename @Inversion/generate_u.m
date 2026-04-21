function obj = generate_u(obj)

% Samples a permeability function

L = obj.cholesky_L_fwd; % Cholesky decomposition
N = length(L);
z = randn(N, 1); % Standard normal random variables
f_sample = L * z; % Sample from the GP (zero-mean)

obj.u_true = f_sample;

end