classdef Volume

    properties
    
        mesh_class;
        darcy_class;

        %% Measure of control volumes in 3D (using thickness approximation)
        volume_measures

        %% Control volume outward facing normals (scaled by edge length)
        volume_outflow_vectors
        

        %% inlet, outlet, and Nuemann boundary node lists
        inlet_nodes;
        outlet_nodes;

        %% connectivity features of the Volume elements
        node_connectivity;
        

        %% legacy code


    end % end properties

    methods
        %% Methods of the Volume class

        function obj = Volume(pressure_class,darcy_class)

            obj.mesh_class = pressure_class.mesh_class;
            obj.darcy_class = darcy_class;

            obj = obj.compute_volume_measures();
            obj = obj.compute_volume_outflow_vectors();


        end

        function obj = compute_volume_measures(obj)
            
            obj.volume_measures = zeros(obj.mesh_class.num_nodes,1);

            for i = 1:obj.mesh_class.num_elements
                
                element_area = obj.mesh_class.element_areas(i);
                control_volumes = obj.mesh_class.elements(i,:);
                sub_volume_measures = element_area/3 * obj.darcy_class.thickness *ones(3,1);
                obj.volume_measures(control_volumes) = obj.volume_measures(control_volumes) + sub_volume_measures;
            end
        end

        function obj = compute_volume_outflow_vectors(obj)

            nodes = obj.mesh_class.nodes;
            elements = obj.mesh_class.elements;
            centroids = obj.mesh_class.centroids;

            %% compute midpoints on triangular elements
            a = (nodes(elements(:,1),:) + nodes(elements(:,2),:))/2;
            b = (nodes(elements(:,2),:) + nodes(elements(:,3),:))/2;
            c = (nodes(elements(:,3),:) + nodes(elements(:,1),:))/2;
            
            %% compute outward triangle element normals
            n1 = [centroids(:,2)-a(:,2) a(:,1)-centroids(:,1)];
            n2 = [centroids(:,2)-b(:,2) b(:,1)-centroids(:,1)];
            n3 = [centroids(:,2)-c(:,2) c(:,1)-centroids(:,1)];
            
            %% outlfow and inflow vectors for each control volume
            n31 = n3 - n1;
            n12 = n1 - n2;
            n23 = n2 - n3;

            %% outward flows for each control volume in this triangle
            obj.volume_outflow_vectors = [n31 n12 n23];

        end


    end % end methods

end