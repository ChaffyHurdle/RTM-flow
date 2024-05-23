classdef CVFEM

    %% Method class to store algorithms, steps and properties of the CVFEM:
    % The cvfem class store the components of the control volume finite
    % element method. This includes classes for the mesh, pressure FEM
    % solver, volume FV solver, darcy flow, visualisation & extra options.

    properties

        mesh_class;
        pressure_class;
        volume_class;
        darcy_class;
        visualise_class;
        options_class;

        time;
        time_step;

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
        inlet_flag;

    end % end properties

    methods
    %% CVFEM class methods:
        % A constructor that takes and stores prebuilt classes of the mesh,
        % pressure, volum, darcy, visualisation and extras classes.

        function obj = CVFEM(mesh_class,pressure_class,volume_class,darcy_class,visualise_class,options_class)

            obj.mesh_class = mesh_class;
            obj.pressure_class = pressure_class;
            obj.volume_class = volume_class;
            obj.darcy_class = darcy_class;
            obj.visualise_class = visualise_class;
            obj.options_class = options_class;

            %% Setting up time and time stepping
            obj.time = 0.0;
            obj.time_step = 0.0;

             %% Setting up active nodes/elements
            obj.inlet_flag = pressure_class.inlet_flag;

            active_nodes = zeros(mesh_class.num_nodes,1);
            active_nodes(obj.inlet_flag==1) = 1;
            
            active_elements = zeros(mesh_class.num_elements,1);
            inlet_idx = find(obj.inlet_flag);
            % At the begining, there are no active elements because no flow moves into
            % the domain through the inlet yet. But in the flux calculation, elements
            % involving inlet nodes needs to be highlighted somehow. 
            % We assigned 0.5 to these elements to distinguish them from finite elements.
            candidate = zeros(mesh_class.num_elements,1);
            for i = 1 : nnz(obj.inlet_flag)
                ival = inlet_idx(i);
                for j = 2 : volume_class.has_node_i(ival,1)+1
                    candidate(volume_class.has_node_i(ival,j))= 1;
                end
            end
            candidate_idx = find(candidate);
            for i = 1 : length(candidate_idx)
                if sum(obj.inlet_flag(mesh_class.elements(candidate_idx(i),:)))==2
                    active_elements(candidate_idx(i)) = 0.5;
                end
            end
            
            obj.inlet_connected_elements = sparse(active_elements==0.5);
            obj.active_nodes = active_nodes;
            obj.active_elements = active_elements;

            %% Setting up volume tracking properties
            obj.volume_filling_times = zeros(mesh_class.num_nodes,1);
            obj.volume_fill_percentage = zeros(mesh_class.num_nodes,1);
            obj.new_filled_volume = [];
            obj = obj.compute_flow_rates();

        end

    end % end methods
end