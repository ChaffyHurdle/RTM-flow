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
        face_normals;

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
                
            obj.nodes = p';
            obj.elements = t(1:4,:)';
            obj.num_nodes = size(p,2);
            obj.num_elements = size(t,2);

            obj = obj.compute_element_faces();
            obj = obj.compute_element_volumes_and_centroids();
            obj = obj.compute_face_areas_and_normals(); 

        end

        function obj = compute_element_faces(obj)

            obj.num_faces = obj.num_elements * 4;

            obj.element_faces = zeros(obj.num_faces,3);

            for i = 1:obj.num_elements

                obj.element_faces(4*(i-1)+1:4*i,:) = nchoosek(obj.elements(i,:),3);
            end

        end

        function obj = compute_element_volumes_and_centroids(obj)

            obj.element_volumes = zeros(obj.num_elements,1);
            obj.centroids = zeros(obj.num_elements,3);

            for i =1:obj.num_elements

                local_nodes = obj.nodes(obj.elements(i,:),:);
                a = local_nodes(1,:); b = local_nodes(2,:); 
                c = local_nodes(3,:); d = local_nodes(4,:);

                product = dot((b-a),cross(c-a,d-a));
                obj.element_volumes(i) = product/6;

                obj.centroids(i,:) = [mean(local_nodes(:,1)) mean(local_nodes(:,2)) mean(local_nodes(:,3))];

            end


        end

        function obj = compute_face_areas_and_normals(obj)

            obj.face_normals = zeros(obj.num_faces,3);
            obj.face_areas = zeros(obj.num_faces,1); 

            for i =1:obj.num_faces

                local_nodes = obj.nodes(obj.element_faces(i,:),:);
                a = local_nodes(1,:); 
                b = local_nodes(2,:); 
                c = local_nodes(3,:);

                normal = cross(c-a,b-a);

                obj.face_areas(i) = 0.5*norm(normal);

                unit_normal = normal./norm(normal);

                if dot(unit_normal,a-obj.centroids(ceil(i/4),:)) < 0

                    unit_normal = -1 * unit_normal;
                 
                end

                obj.face_normals(i,:) = unit_normal;

            end

        end

    end

end