classdef Pressure

    properties
    
        mesh_class;
        inlet_func;
        vent_func;
        p_D;

        %% FEM system of equations
        stiffness_matrix;
        load_vector;

        %% FEM solution of pressure problem
        pressure;
        pressure_gradient;

        %% inlet, outlet, and Nuemann boundary node lists
        inlet_nodes;
        outlet_nodes;
        Neumann_nodes;
        free_nodes;

    end % end properties

    methods

        function obj = Pressure(mesh_class,inlet_script,outlet_script,p_D)

            obj.mesh_class = mesh_class;
            obj.inlet_func = str2func(inlet_script);
            obj.vent_func = str2func(outlet_script);
            obj.p_D = p_D;

            num_dofs = mesh_class.num_nodes;
            obj.pressure = obj.p_D(mesh_class.nodes,[],[]);
            obj.pressure_gradient = zeros(mesh_class.num_elements,2);

            obj.stiffness_matrix = spalloc(num_dofs,num_dofs,10*num_dofs);
            obj.load_vector = zeros(num_dofs,1);

        end

    end % end methods

end