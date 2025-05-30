%% Set sqrtN^2 (i.e. N) sensor locs for square grid (CASE 1, CASE 2)
sqrtN = 3;
sensor_locs_x = 1/(2*sqrtN) + linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = sensor_locs_x;
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];
save(strcat("Case1/meshes_sensors/sensor_locs_",num2str(sqrtN^2),".mat"),"sensor_locs")
save(strcat("Case2/meshes_sensors/sensor_locs_",num2str(sqrtN^2),".mat"),"sensor_locs")

%% Set 2*sqrtN^2 (i.e. 2N) sensor locs for fork (CASE 3)
sqrtN = 10;
sensor_locs_x = 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y = 0.25 + 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
[sensor_locs_x,sensor_locs_y] = meshgrid(sensor_locs_x,sensor_locs_y);
sensor_locs_x = reshape(sensor_locs_x,[],1);
sensor_locs_y = reshape(sensor_locs_y,[],1);
sensor_locs = [sensor_locs_x sensor_locs_y];

sensor_locs_x_prong = 0.5 + 0.5/(2*sqrtN) + 0.5*linspace(0,sqrtN-1,sqrtN)/sqrtN;
sensor_locs_y_prong = sort(unique(sensor_locs_y));
for i = 1:sqrtN
    for j = 1:sqrtN/2
        %prong_upper_y = 0.5*sensor_locs_x_prong(i) + (sensor_locs_y_prong(j)-0.25);
        prong_lower_y = -0.5*sensor_locs_x_prong(i) + (sensor_locs_y_prong(j)+0.25);
        sensor_locs = [sensor_locs; sensor_locs_x_prong(i) 1-prong_lower_y];
        sensor_locs = [sensor_locs; sensor_locs_x_prong(i) prong_lower_y];
    end
end
save(strcat("Case3/meshes_sensors/sensor_locs_",num2str(2*sqrtN^2),".mat"),"sensor_locs")


%% Set n_per_turn*n_radii sensor locs for torus (CASE 4)
n_per_turn = 6;
n_radii = 4;
turns = 2*pi * (0:(n_per_turn-1)) / n_per_turn;
radii = 0.25 + 0.75/(2*n_radii) + 0.75*linspace(0,n_radii-1,n_radii)/n_radii;
sensor_locs_radii = reshape(repmat(radii,n_per_turn,1),1,[]);
sensor_locs_theta = repmat(turns,1,n_radii) + reshape(repmat(double(~mod(1:n_radii,2))*(turns(2)-turns(1))/2,n_per_turn,1),1,[]);
sensor_locs_x = sensor_locs_radii .* cos(sensor_locs_theta);
sensor_locs_y = sensor_locs_radii .* sin(sensor_locs_theta);
sensor_locs = [sensor_locs_x' sensor_locs_y'];
save(strcat("Case4/meshes_sensors/sensor_locs_",num2str(n_per_turn*n_radii),".mat"),"sensor_locs")