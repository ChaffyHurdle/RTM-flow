function Q_i = compute_representer_ij(grad_lambdas,obj_data)

grad_pressures = obj_data.RTMflow_class.pressure_gradients;
active_elements = obj_data.RTMflow_class.all_active_elements;
times = obj_data.RTMflow_class.times;
diff_times = diff(times);
Q_i = zeros(1,obj_data.mesh_class.num_elements);
u = obj_data.u;

    for i = 1:(length(times)-1)
        f_i = active_elements(:,i) .* sum(grad_lambdas(:,:,i) .* grad_pressures{i},2);
        f_iplus1 = active_elements(:,i+1) .* sum(grad_lambdas(:,:,i+1) .* grad_pressures{i+1},2);
        Q_i = Q_i + diff_times(i) * (f_i + f_iplus1)'/2;
    end
    Q_i = - Q_i .* exp(u);

    Q_i(abs(Q_i)<1e-10)=0;
end