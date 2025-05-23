function [boldR,curlyR,c] = compute_representers(u,ups,C,params,t)

M = params.nsensors;
N = params.nobtimes;
sensor_locs = params.sensor_locs;

Q_mat = zeros(params.Nx,M*t);

switched_on_inds = zeros(1,M*t);
steady_state_inds = zeros(1,M*t);

for k = 1:M*t
    i = mod(k-1,M)+1;
    j = ceil(k/M);
    x_i = sensor_locs(i);
    Ups_j = ups(j);

    switched_on = (x_i <= Ups_j);
    switched_on_inds(k) = switched_on;

    steady_state = (Ups_j == params.L);
    steady_state_inds(k) = steady_state;

    if switched_on
        premultiplier = (params.p_I - params.p_0)/F_u(Ups_j,u,params)^2;
        [~,nearest_x_to_ups] = min(abs(params.x_locations - Ups_j));
        first = F_u(Ups_j,u,params)*exp(-u).*(params.x_locations < x_i);
        second = F_u(x_i,u,params)*exp(-u).*(params.x_locations < Ups_j);
        third_1 = (F_u(x_i,u,params)*exp(-u(nearest_x_to_ups)))/F_u(Ups_j,u,params);
        third_2 = exp(-u) .* max(Ups_j - params.x_locations,0);
        
        Q_i = premultiplier * (first - second + switched_on*third_1*third_2);
        Q_mat(:,k) = Q_i;
    end
end
boldR = params.dx * C * Q_mat;

curlyR = zeros(t*M,t*M);
c = zeros(t*M,1);
for i = 1:t*M
    for j = i:t*M
        curlyR(i,j) = params.dx * sum(Q_mat(:,i) .* boldR(:,j));
        curlyR(j,i) = curlyR(i,j);
    end
    c(i) = params.dx * sum(Q_mat(:,i) .* (params.u_bar - u)');
end

end