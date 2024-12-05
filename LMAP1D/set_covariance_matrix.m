function C_mat = set_covariance_matrix(params)

C_mat = zeros(params.Nx,params.Nx);
nu = params.nu;
l = params.l;
var_matern = params.var_matern;

for i = 1:params.Nx
    for j = i:params.Nx
        if i ~= j
            h = abs(params.x_locations(i)-params.x_locations(j));
            part1 = 2^(1-nu)/gamma(nu);
            part2 = (sqrt(2*nu)*h/l)^nu;
            part3 = besselk(nu,sqrt(2*nu)*h/l);
            c_ij = var_matern * part1 * part2 * part3;
            C_mat(i,j) = c_ij;
            C_mat(j,i) = c_ij;
        else
            C_mat(i,j) = var_matern;
        end
    end
end

end