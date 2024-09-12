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
        p_tilde;
        grad_p_tilde;
        is_moving_boundary;
        dirichlet_nodes_matrix;
        bndry_cond;
        v_h;
        bndry_conds;
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
            obj.bndry_cond = zeros(length(mesh_class.nodes),1);
            obj.v_h = zeros(length(mesh_class.nodes),1);
            obj.bndry_conds = [];
        end
    end
end