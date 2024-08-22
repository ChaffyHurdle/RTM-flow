function h = generate_h(obj)

L = obj.cholesky_L_fwd; % Cholesky decomposition
NL = length(L);
z = randn(NL, 1); % Standard normal random variables
f_sample = L * z; % Sample from the GP (zero-mean)

h = f_sample;

end