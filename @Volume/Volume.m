classdef Volume

    properties
    
        mesh_class;

        %% FEM system of equations
        

        %% FEM solution of pressure problem
        

        %% inlet, outlet, and Nuemann boundary node lists
        inlet_nodes;
        outlet_nodes;

    end % end properties

    methods

        function obj = Volume(pressure_class)

            obj.mesh_class = pressure_class.mesh_class;


        end

    end % end methods

end