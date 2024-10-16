function J = evaluate_cost_function(obj,u,RTM_class)

C0_minushalf = obj.inverse_class.C0_minushalf;
Sigma = reshape(obj.inverse_class.Sigma,[],1);
Sigma_minus_half = diag(1./sqrt(Sigma));
data = reshape(obj.inverse_class.data,[],1);
Gu = reshape(RTM_class.pressure_data,[],1);

%prior_contribution_2 = sum((C0_minushalf*(u - obj.u0)').^2 .* obj.mesh_class.element_areas)
prior_contribution = norm(C0_minushalf*(u - obj.u0)')^2;
lhood_contribution = norm(Sigma_minus_half*(data - Gu))^2;

J = prior_contribution + lhood_contribution;

end