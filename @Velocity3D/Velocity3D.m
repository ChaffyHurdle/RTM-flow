classdef Velocity3D

    properties
    
        Delaunay_mesh_class;
        physics_class;

        %% piecewise constant velocity field
        time;
        velocity;

    end % end properties

    methods
        %% Methods of the Volume class

        function obj = Velocity3D(Delaunay_mesh_class,physics_class)

            obj.Delaunay_mesh_class = Delaunay_mesh_class;
            obj.physics_class = physics_class;
            
            num_elem = Delaunay_mesh_class.num_elements;
            dim = size(Delaunay_mesh_class.nodes,2);

            obj.time = 0;
            obj.velocity = zeros(num_elem,dim);
        end    

    end % end methods

end