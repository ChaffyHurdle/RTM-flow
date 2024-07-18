classdef Pressure

    properties
    
        %% mesh and physics classes stored for easy access
        mesh_class;
        physics_class;

        %% functions to determine inlet/vent pressures
        inlet_func;
        vent_func;
        p_D;

        %% time stored for reference
        time;

        %% FEM system of equations
        stiffness_matrix;
        load_vector;

        %% FEM solution of pressure problem
        pressure;
        shape_fun_gradients;
        pressure_gradient;
        
        %% inlet, outlet, and Nuemann boundary node lists
        is_inlet;
        is_Dirichlet;
        is_vent;
        is_Neumann;

        %% degrees of freedom in the FEM system
        is_node_active;
        new_active_elements;


    end % end properties

    methods

        function obj = Pressure(mesh_class,physics_class)

            obj.mesh_class = mesh_class;
            obj.physics_class = physics_class;

            % Hard code inlet/vent positions
            %obj.inlet_func = inlet_func;
            obj.inlet_func = @(x) (x(1) == 0);
            %obj.vent_func = vent_func;
            obj.vent_func = @(x) (x(1) == 1);
            %obj.p_D = p_D;

            %% set time to zero
            obj.time = 0.0;

            %% Computing information for boundary conditions
            obj = obj.compute_inlets_outlets();

            %% Allocating pressure and pressure gradient
            num_dofs = mesh_class.num_nodes;

            % Hard code Dirichlet condition
            %obj.pressure = obj.pressure.p_D(obj);
            obj.pressure = physics_class.p_0*ones(num_dofs,1);
            obj.pressure(obj.is_inlet) = physics_class.p_I;

            obj = obj.compute_shape_fun_gradients();
            obj = obj.compute_pressure_gradient();

            %% Allocating FEM system of equations
            obj.stiffness_matrix = spalloc(num_dofs,num_dofs,10*num_dofs);
            obj.load_vector = zeros(num_dofs,1);

        end

    end % end methods

end