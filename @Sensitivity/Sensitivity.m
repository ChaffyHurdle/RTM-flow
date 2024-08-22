classdef Sensitivity

    properties

        % Classes
        mesh_class;
        darcy_class_u;
        
        % log-permeability and perturbation
        u;
        h;

        % Forward solves
        RTMflow_u;
        RTMflow_u_plus_h;
    end

    methods
        function obj = Sensitivity(mesh_class,darcy_class_u,u,h)
            
            % Classes
            obj.mesh_class= mesh_class;
            obj.darcy_class_u = darcy_class_u;

            % log-permeabilities
            obj.u = u;
            % perturbation
            obj.h = h;
        end
    end
end