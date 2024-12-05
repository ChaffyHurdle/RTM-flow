function [d,Sigma] = set_data(p,sigma1,sigma2,params)

M = params.nsensors;
N = params.nobtimes;
Sigma1 = ( sigma1 * abs(p) ) .^ 2;
Sigma2 = ( sigma2 * abs(max(p,[],"all")-min(p,[],"all")) ) .^ 2;
Sigma = Sigma1 + Sigma2;
noise = sqrt(Sigma1) .* normrnd(0,1,M,N) + sqrt(Sigma2) .* normrnd(0,1,M,N);
d = p + noise;

end