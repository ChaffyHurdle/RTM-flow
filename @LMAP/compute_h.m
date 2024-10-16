function h = compute_h(obj)

prior_dist = obj.u0 - obj.u;
linearised_data_misfit = reshape(obj.inverse_class.data,[],1) ...
                       - reshape(obj.RTMflow_class.pressure_data,[],1) - obj.d/(1+obj.alpha);
RplusSigma = obj.tildePmat + (1+obj.alpha)*diag(reshape(obj.inverse_class.Sigma,[],1));
h = prior_dist/(1+obj.alpha) + (obj.R*(RplusSigma\linearised_data_misfit))';

end