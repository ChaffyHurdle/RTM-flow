function obj = construct_linear_system_t(obj)

obj.tildePmat = obj.Q' * (obj.R .* obj.mesh_class.element_areas);
obj.d = obj.Q' * ((obj.u0 - obj.u)'.*obj.mesh_class.element_areas);

end