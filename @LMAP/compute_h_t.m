function h = compute_h_t(obj,t,n_fwds)

% Computes n_fwds proposals for h with [\alpha,\alpha*k,\alpha*k^2,\alpha*k^3,\alpha*k^4]
h = zeros(obj.mesh_class.num_elements,n_fwds);
prior_dist = obj.u0 - obj.u;
curr_alpha = obj.alpha;

for i = 1:n_fwds
    alpha = curr_alpha*obj.scale^(i-1);
    linearised_data_misfit = reshape(obj.inverse_class.data(:,1:t),[],1) ...
                           - reshape(obj.RTMflow_class.pressure_data(:,1:t),[],1) - obj.d/(1+alpha);
    RplusSigma = obj.tildePmat + (1+alpha)*diag(reshape(obj.inverse_class.Sigma(:,1:t),[],1));
    h(:,i) = prior_dist/(1+alpha) + (obj.R*(RplusSigma\linearised_data_misfit))';
end

end