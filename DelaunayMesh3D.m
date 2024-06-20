classdef DelaunayMesh3D

    %% properties
    properties
        %% Mesh Data
        nodes;
        elements;
        boundary_nodes;
        boundary_faces;
        centroids;

        %% Element properties
        element_volumes;
        element_faces;
        face_areas;

        %% Mesh Counting Properties
        num_nodes;
        num_elements;
        num_faces;
        num_boundary_nodes;
        num_boundary_faces;
    end

    %% methods
    methods

        function obj = DelaunayMesh3D(gmsh_filename,h_max)

            %% extract gmsh mesh properties
            run(gmsh_filename); %runs matlab script output is msh
            p = msh.POS'; 
            t = msh.TETS';

            if nargin == 2

                %% matlab computes new mesh
                model = createpde;
                model.geometryFromMesh(p,t(1:4,:));

                mesh = model.generateMesh(Hmax=h_max,GeometricOrder="linear");
                [p,e,t] = meshToPet(mesh);
            end
                
            obj.nodes = p;
            obj.elements = t(1:4,:);
            obj.num_nodes = size(p,2);
            obj.num_elements = size(t,2);

            obj = obj.compute_element_faces();
            obj = obj.compute_element_volumes();

        end

        function obj = compute_element_faces(obj)

            num_faces = obj.num_elements * 4;

            obj.element_faces = zeros(3,num_faces);

            for i = 1:obj.num_elements

                obj.element_faces{i} = nchoosek(obj.elements(:,i),3);
            end


            obj.num_faces = num_faces;
        end

        function obj = compute_element_volumes(obj)

            obj.element_volumes = zeros(obj.num_elements,1);

            for i =1:obj.num_elements

                local_nodes = obj.nodes(:,obj.elements(:,i));
                a = local_nodes(:,1); b = local_nodes(:,2); 
                c = local_nodes(:,3); d = local_nodes(:,4);

                product = abs((b-a)' * cross(c-a,d-a));
                obj.element_volumes(i) = product/6;

            end


        end

        function obj = compute_face_areas()
            


        end

    end

end