function h = compute_h_t(obj,curly_R,bold_R,d,t)

h = zeros(obj.mesh_class.num_elements,5);
prior_dist = obj.u0 - obj.u;
curr_alpha = obj.alpha;

for i = 1:5
    alpha = curr_alpha/obj.scale^(i-1);
    linearised_data_misfit = reshape(obj.inverse_class.data(:,1:t),[],1) ...
                           - reshape(obj.RTMflow_class.pressure_data(:,1:t),[],1) - d/(1+alpha);
    RplusSigma = curly_R + (1+alpha)*diag(reshape(obj.inverse_class.Sigma(:,1:t),[],1));
    
    h(:,i) = prior_dist/(1+alpha) + (bold_R*(RplusSigma\linearised_data_misfit))';
end

end