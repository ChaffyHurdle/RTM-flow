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
        is_moving_boundary; 
        is_moving_boundary_u_plus_h;
        is_moving_boundary_ob_times;
        is_moving_boundary_ob_times_u_plus_h;
        dirichlet_nodes_matrix;
        pressures_u;
        pressures_u_plus_h;
        active_nodes_u;
        time_inds_u;
        time_inds_u_plus_h;

        % Sensitivity objects
        p_tilde; % p_tilde at every time
        p_tildes; % p_tilde at observation times
        p_tildes_at_sensors; % p_tilde at observation times and sensor locs
        grad_p_tilde; % gradient of p_tilde at time snapshot
        bndry_cond; % value of p_tilde on boundary
        v_h; % normal perturbation at boundary
        bndry_conds; % value of p_tilde on boundary at all times
        edge_data; % data describing moving front position at fixed time
        all_edge_data; % data describing moving front position at all times
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
            obj.all_edge_data = cell(1);
        end
    end
end