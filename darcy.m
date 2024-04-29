classdef darcy
    
    %% Darcy class to Darcy's law and important material properties:
    % The darcy class serves to store the domain material properties:
    % viscosity, porosity, thickness (2d -> 3d scaling). Each variable is
    % a matlab function which inputs a point coordinate and time.
    % 
    % viscosity is ....
    % porosity is ....
    % thickness is ...

    properties
        %% Material properties
        viscosity;
        porosity;
        thickness;
    end

    methods

        %% Darcy class methods:
        % A constructor that inputs 3 functions for each of the material
        % properties. These must be defined within matlab before entering
        % into the constructor function.

        function obj = darcy(viscosity, porosity, thickness)

            obj.viscosity = viscosity;
            obj.porosity = porosity;
            obj.thickness = thickness;

        end


    end


end
