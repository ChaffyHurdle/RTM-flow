classdef Volume

    properties
    
        mesh_class;
        darcy_class;

        %% Measure of control volumes in 3D (using thickness approximation)
        volume_measures;

        %% Control volume outward facing normals(multiplied by edge length)
        volume_outflow_vectors;

        %% Inlet, outlet, and Nuemann boundary node lists
        inlet_nodes;
        outlet_nodes;

        %% Connectivity features of the Volume elements
        node_connectivity;
        element_connectivity;

        %% Tracking volumes filled
        volume_fill_percentage;
        volume_filling_times;
        volume_rates_of_flow;
        new_filled_volume;

        %% Features for the FEM pressure solver
        inlet_connected_elements;
        active_nodes;
        active_elements;

        %% Legacy code
        has_node_i;
        bndry_nodes;
        nb_nodes;

        inlet_flag;

    end % end properties

    methods
        %% Methods of the Volume class

        function obj = Volume(pressure_class,darcy_class)

            obj.mesh_class = pressure_class.mesh_class;
            obj.darcy_class = darcy_class;

            %% Set initial fill times and percentages
            obj.volume_filling_times = zeros(obj.mesh_class.num_elements,1);
            obj.volume_fill_percentage = zeros(obj.mesh_class.num_elements,1);

            obj = obj.compute_volume_measures();
            obj = obj.compute_volume_outflow_vectors();
            obj = obj.compute_connectivity();

        end    

    end % end methods

end