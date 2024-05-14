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

            obj.time = 0.0;

        end

    end % end methods
end