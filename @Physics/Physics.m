classdef Physics
    
    %% Physics class to store Darcy's law and important material properties:
    % The physics class serves to store the domain material properties:
    % viscosity, porosity, thickness (2d -> 3d scaling).

    properties
        %% Material properties
        viscosity;
        porosity;
        thickness;
        p_I;
        p_0;
        permeability;
        sensor_locs;
        nsensors;
        observation_times;
        nobservations;
        T;
    end

    methods

        %% Darcy class methods:
        % A constructor that inputs 3 functions for each of the material
        % properties. These must be defined within matlab before entering
        % into the constructor function.

        function obj = Physics(viscosity, porosity, thickness, inlet_pressure, outlet_pressure, permeability, sensor_locs, observation_times, T)

            obj.viscosity = viscosity;
            obj.porosity = porosity;
            obj.thickness = thickness;
            obj.p_I = inlet_pressure;
            obj.p_0 = outlet_pressure;
            obj.permeability = permeability;
            obj.sensor_locs = sensor_locs;
            obj.observation_times = observation_times;
            obj.T = T;
            obj.nsensors = length(sensor_locs);
            obj.nobservations = length(observation_times);

        end


    end

end
