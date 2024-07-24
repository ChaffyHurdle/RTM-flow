function obj = generate_u(obj)

L = obj.cholesky_L_fwd; % Cholesky decomposition
NL = length(L);
z = randn(NL, 1); % Standard normal random variables
f_sample = L * z; % Sample from the GP (zero-mean)

obj.u_true = f_sample;

end