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
        p_tildes;
        p_tildes_at_sensors;
        grad_p_tilde;
        is_moving_boundary;
        is_moving_boundary_ob_times;
        is_moving_boundary_u_plus_h;
        is_moving_boundary_ob_times_u_plus_h
        dirichlet_nodes_matrix;
        bndry_cond;
        v_h;
        bndry_conds;
        pressures_u;
        pressures_u_plus_h;
        active_nodes_u;
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
            obj.v_h = zeros(mesh_class.num_elements,2);
            obj.bndry_conds = [];
            obj.grad_p_tilde = zeros(mesh_class.num_elements,2);
        end
    end
end