function obj = generate_data(obj,pressure_data,sigma0)

obj.pressure_noiseless = pressure_data;

pressure_transformed = pressure_data - min(pressure_data);

obj.Sigma = sigma0^2*(max(pressure_transformed,[],'all') ...
    - min(pressure_transformed,[],'all'))^2*ones(size(pressure_data));
noise = normrnd(0,1,size(pressure_data)).*sqrt(obj.Sigma);
obj.data = pressure_data + noise;

end