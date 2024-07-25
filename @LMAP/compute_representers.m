function obj = compute_representers(obj)
disp("Computing representers")

% Various init variables
nobs = obj.physics_class.nobservations;
nsensors = obj.physics_class.nsensors;
representers = zeros(obj.mesh_class.num_elements,nsensors*nobs);
times = obj.RTMflow_class.times;
diff_times = diff(times);
disp("Moving grad lambdas to pool")
grad_lambdas = parallel.pool.Constant(obj.grad_lambdas);
disp("Moving grad pressures to pool")
grad_pressures = parallel.pool.Constant(obj.RTMflow_class.pressure_gradients);
permeability = obj.physics_class.permeability;
disp("Moving active_elements to pool")
active_elements = parallel.pool.Constant(obj.RTMflow_class.all_active_elements);


parfor k = 1:nobs * nsensors
    grad_lambda_i = grad_lambdas.Value(:,:,:,k);
    grad_pressures_i = grad_pressures.Value;
    active_elements_i = active_elements.Value;
    disp(k)

    R_i = compute_representer_i(times,diff_times, ...
        grad_lambda_i,grad_pressures_i,active_elements_i,permeability);
    representers(:,k) = R_i;
end
obj.representers = representers;

end


function R_i = compute_representer_i(times,diff_times, ...
        grad_lambda_i,grad_pressures_i,active_elements_i,permeability)
    R_i = zeros(1,length(grad_lambda_i));

    for i = 1:(length(times)-1)
        f_i = active_elements_i(:,i) .* sum(grad_lambda_i(:,:,i) .* grad_pressures_i{i},2);
        f_iplus1 = active_elements_i(:,i+1) .* sum(grad_lambda_i(:,:,i+1) .* grad_pressures_i{i+1},2);
        R_i = R_i + diff_times(i) * (f_i + f_iplus1)'/2;
    end
    R_i = - R_i .* permeability;
end