function [J,lhood_contribution] = evaluate_cost_function_t(obj,u,RTM_class,t)

C0_minushalf = obj.inverse_class.C0_minushalf;
Sigma = reshape(obj.inverse_class.Sigma(:,1:t),[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
data = reshape(obj.inverse_class.data(:,1:t),[],1);
Gu = reshape(RTM_class.pressure_data(:,1:t),[],1);

prior_contribution = norm(C0_minushalf*(u - obj.u0)')^2;
lhood_contribution = norm(Sigma_minus_half*(data - Gu))^2;

J = 0.5*(prior_contribution + lhood_contribution);

end