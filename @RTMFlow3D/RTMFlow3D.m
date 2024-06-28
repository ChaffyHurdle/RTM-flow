classdef RTMFlow3D

    %% 3D version of RTM flow

    properties

        Delaunay_mesh_class;
        pressure_class;
        velocity_class;
        physics_class;
        visualise_class;

        time;
        time_step;

        %% Tracking volumes filled
        volume_fill_percentage;
        volume_filling_times;
        volume_rates_of_flow;
        new_filled_volume;
        active_elements;

    end % end properties

    methods
    %% CVFEM class methods

    function obj = RTMFlow3D(Delaunay_mesh_class,physics_class,pressure_class)
            
            %% store variables
            obj.Delaunay_mesh_class = Delaunay_mesh_class;
            obj.pressure_class = pressure_class;
            obj.physics_class = physics_class;
            obj.velocity_class = Velocity3D(obj.Delaunay_mesh_class,physics_class);

            %% Default to no visuals
            obj.visualise_class = Visualisation();

            %% Setting up time and time stepping
            obj.time = 0.0;
            obj.time_step = 0.0;

            %% Setting up active nodes/elements
            inlet_flag = pressure_class.is_inlet;
            obj.active_elements = zeros(Delaunay_mesh_class.num_elements,1);
            inlet_nodes = find(inlet_flag);

            is_inlet_connected = ismember(obj.Delaunay_mesh_class.elements, inlet_nodes);
            active_element_list = find(sum(is_inlet_connected,2)>0);
            obj.active_elements(active_element_list) = 0.5;

            %% Setting up volume tracking properties
            obj.volume_filling_times = zeros(Delaunay_mesh_class.num_elements,1);
            obj.volume_fill_percentage = zeros(Delaunay_mesh_class.num_elements,1);
            obj.new_filled_volume = [];
            obj = obj.compute_flow_rates();

        end

    end % end methods
end