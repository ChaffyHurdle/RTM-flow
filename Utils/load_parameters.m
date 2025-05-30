function [matern_args,mu,phi,thickness,p_I,p_0,observation_times,T] = load_parameters(case_study)

var_matern = 0.25;
nu_matern = 1.5;

switch case_study
    case 1
        length_scale = 0.1;
        mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
    case 2
        length_scale = 0.1;
        mu = 1; phi = 0.5; thickness = 1; p_I = 6; p_0 = 1;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
    case 3
        length_scale = 0.1;
        mu = 0.5; phi = 1; thickness = 1; p_I = 1; p_0 = 0;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
    case 4
        length_scale = 0.2;
        mu = 0.5; phi = 1; thickness = 1; p_I = 2; p_0 = 1;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
end

matern_args = [var_matern,length_scale,nu_matern];


end