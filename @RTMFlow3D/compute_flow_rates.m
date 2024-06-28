function obj = compute_flow_rates(obj)

%% Extract needed information
fFactor = obj.volume_fill_percentage;
obj.volume_rates_of_flow = zeros(obj.Delaunay_mesh_class.num_elements,1);

%% construct a list of possible elements to be added
candidate_elem = find(~xor(obj.active_elements>=0.5,1-fFactor > 1e-15));
num_candidates = length(candidate_elem);
is_candidate = false(num_candidates,1);

for i = 1:num_candidates

    elem_num = candidate_elem(i);
    face_nums = (4*(elem_num-1)+1:4*elem_num)';
    other_faces = obj.Delaunay_mesh_class.face_connectivity(face_nums);
    
    %% if candidate is face connected to an active element, flag for flow
    for j = 1:4

        other_elem_num = ceil(other_faces(j)/4);

        if other_elem_num == 0
            continue
        end

        if obj.active_elements(other_elem_num) >= 0.5
            is_candidate(i) = true;
            break
        end
    end

end

final_candidates = candidate_elem(is_candidate);

for i = 1:length(final_candidates)
    
    elem_num = final_candidates(i);
    face_nums = (4*(elem_num-1)+1:4*elem_num)';

    local_velocity = obj.velocity_class.velocity(elem_num,:);
    
    unit_normals = obj.Delaunay_mesh_class.face_normals(face_nums,:);
    connected_faces = obj.Delaunay_mesh_class.face_connectivity(face_nums);

    flux = 0;
    
    for j = 1:4

        if connected_faces(j) == 0
            other_velocity = [0 0 0]; % if boundary no flux contribution
        else
            other_velocity = obj.velocity_class.velocity(ceil(connected_faces(j)/4),:);
        end
        
        flux = flux + dot(local_velocity-other_velocity,unit_normals(j,:));
    end

    obj.volume_rates_of_flow(elem_num) = flux;

end

end