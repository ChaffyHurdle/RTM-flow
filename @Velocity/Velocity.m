classdef Velocity

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

        function obj = Velocity(pressure_class,darcy_class)

            obj.mesh_class = pressure_class.mesh_class;
            obj.darcy_class = darcy_class;

            %% Set initial fill times and percentages
            obj.volume_filling_times = zeros(obj.mesh_class.num_elements,1);
            obj.volume_fill_percentage = zeros(obj.mesh_class.num_elements,1);

            obj = obj.compute_volume_measures();
            obj = obj.compute_volume_outflow_vectors();
            obj = obj.compute_connectivity();

            %% Setting up active nodes/elements
            obj.inlet_flag = pressure_class.inlet_flag;

            active_nodes = zeros(obj.mesh_class.num_nodes,1);
            active_nodes(obj.inlet_flag==1) = 1;
            
            active_elements = zeros(obj.mesh_class.num_elements,1);
            inlet_idx = find(obj.inlet_flag);
            % At the begining, there are no active elements because no flow moves into
            % the domain through the inlet yet. But in the flux calculation, elements
            % involving inlet nodes needs to be highlighted somehow. 
            % We assigned 0.5 to these elements to distinguish them from finite elements.
            candidate = zeros(obj.mesh_class.num_elements,1);
            for i = 1 : nnz(obj.inlet_flag)
                ival = inlet_idx(i);
                for j = 2 : obj.has_node_i(ival,1)+1
                    candidate(obj.has_node_i(ival,j))= 1;
                end
            end
            candidate_idx = find(candidate);
            for i = 1 : length(candidate_idx)
                if sum(obj.inlet_flag(obj.mesh_class.elements(candidate_idx(i),:)))==2
                    active_elements(candidate_idx(i)) = 0.5;
                end
            end
            
            obj.inlet_connected_elements = sparse(active_elements==0.5);
            obj.active_nodes = active_nodes;
            obj.active_elements = active_elements;

        end    

    end % end methods

end