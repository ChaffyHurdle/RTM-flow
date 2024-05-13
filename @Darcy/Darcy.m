classdef Darcy
    
    %% Darcy class to Darcy's law and important material properties:
    % The darcy class serves to store the domain material properties:
    % viscosity, porosity, thickness (2d -> 3d scaling). Each variable is
    % a matlab function which inputs an arg structure containing the 
    % necessary information.
    % 
    % viscosity is ....
    % porosity is ...
    % thickness is ...
    % permeability is ...

    properties
        %% Material properties
        viscosity;
        porosity;
        thickness;
        permeability;
    end

    methods

        %% Darcy class methods:
        % A constructor that inputs 3 functions for each of the material
        % properties. These must be defined within matlab before entering
        % into the constructor function.

        function obj = Darcy(viscosity, porosity, thickness, permeability)

            obj.viscosity = viscosity;
            obj.porosity = porosity;
            obj.thickness = thickness;
            obj.permeability = permeability;

        end


    end

end
