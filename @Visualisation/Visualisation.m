classdef Visualisation

    properties

        %% booleans to decide what to plot
        is_plotting_pressure;
        is_plotting_velocity;
        is_plotting_volume_flow;
        is_plotting_volume;

        %% booleans to decide what to annimate
        is_animate_pressure;
        is_animate_velocity;
        is_animate_volume_flow;
        is_animate_volume;

    end % end properties

    methods

        function obj = Visualisation()
            
            %% default is only plot the volume
            obj.is_plotting_pressure = true;
            obj.is_plotting_velocity = false;
            obj.is_plotting_volume = true;
    
            %% no animations by default (heavy on cpu)
            obj.is_animate_pressure = false;
            obj.is_animate_velocity = false;
            obj.is_animate_volume_flow = false;
            obj.is_animate_volume = false;

        end

     
    end
end