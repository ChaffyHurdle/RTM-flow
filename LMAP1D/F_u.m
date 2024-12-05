function Fux = F_u(x,u,params)

H = 0.5 + 0.5 * tanh(1000*(x(:)-params.x_locations(:)'));
Fux = params.dx * sum((H * (exp(-u))'),2);

end
    