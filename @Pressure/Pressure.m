classdef Pressure

    properties
    
        mesh_class;
        darcy_class;
        inlet_func;
        vent_func;
        p_D;

        %% FEM system of equations
        stiffness_matrix;
        load_vector;

        %% FEM solution of pressure problem
        pressure;
        shape_fun_gradients;
        pressure_gradient;
        

        %% inlet, outlet, and Nuemann boundary node lists
        inlet_nodes;
        outlet_nodes;
        Neumann_nodes;
        free_nodes;

        %% Legacy code
        bndry_nodes;
        nb_nodes;

        inlet_flag;
        inlet_pos;
        Dirichlet;
        vent_flag;
        vent_idx;
        vent_elem;
        Neumann_flag;

    end % end properties

    methods

        function obj = Pressure(mesh_class,darcy_class,inlet_func,vent_func,p_D)

            obj.mesh_class = mesh_class;
            obj.darcy_class = darcy_class;
            obj.inlet_func = inlet_func;
            obj.vent_func = vent_func;
            obj.p_D = p_D;

            %% Allocating pressure and pressure gradient
            obj.pressure = obj.p_D(mesh_class.nodes,[],[],[]);
            obj = obj.compute_shape_fun_gradients();
            obj = obj.compute_pressure_gradient();

            %% Allocating FEM system of equations
            num_dofs = mesh_class.num_nodes;
            obj.stiffness_matrix = spalloc(num_dofs,num_dofs,10*num_dofs);
            obj.load_vector = zeros(num_dofs,1);

            %% Computing information for boundary conditions
            obj = obj.compute_inlets_outlets();
            

        end

    end % end methods

end