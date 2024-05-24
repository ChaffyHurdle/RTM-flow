classdef Volume

    properties
    
        mesh_class;

        %% Measure of control volumes in 3D (using thickness approximation)
        volume_measures;

        %% Control volume outward facing normals(multiplied by edge length)
        volume_outflow_vectors;

        %% Connectivity and boundary features of the Volume elements
        node_connectivity;
        element_connectivity;
        has_node_i;
        bndry_nodes;
        nb_nodes;

    end % end properties

    methods
        %% constructor of the Volume class

        function obj = Volume(mesh_class,darcy_class)

            obj.mesh_class = mesh_class;

            obj = obj.compute_volume_measures(darcy_class.thickness);
            obj = obj.compute_volume_outflow_vectors();
            obj = obj.compute_connectivity();

        end    

    end % end methods

end