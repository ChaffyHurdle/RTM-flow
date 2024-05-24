classdef Pressure

    properties
    
        %% mesh and darcy classes stored for easy access
        mesh_class;
        darcy_class;

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
        inlet_flag;
        inlet_pos;
        Dirichlet;
        vent_flag;
        vent_idx;
        vent_elem;
        Neumann_flag;

        %% degrees of freedom in the FEM system
        active_nodes;
        new_active_elements;


    end % end properties

    methods

        function obj = Pressure(mesh_class,darcy_class,inlet_func,vent_func,p_D)

            obj.mesh_class = mesh_class;
            obj.darcy_class = darcy_class;
            obj.inlet_func = inlet_func;
            obj.vent_func = vent_func;
            obj.p_D = p_D;

            %% set time to zero
            obj.time = 0.0;

            %% Allocating pressure and pressure gradient
            num_dofs = mesh_class.num_nodes;
            obj.pressure = zeros(num_dofs,1);
            obj = obj.compute_shape_fun_gradients();
            obj = obj.compute_pressure_gradient();

            %% Allocating FEM system of equations
            obj.stiffness_matrix = spalloc(num_dofs,num_dofs,10*num_dofs);
            obj.load_vector = zeros(num_dofs,1);

            %% Computing information for boundary conditions
            obj = obj.compute_inlets_outlets();
            

        end

    end % end methods

end