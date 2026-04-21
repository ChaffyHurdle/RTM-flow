function [matern_args,mu,phi,thickness,p_I,p_0,observation_times,T] = load_parameters(case_study)

% Parameters for each case
switch case_study
    case 1 % square
        length_scale = 0.1; var_matern = 0.25; nu_matern = 1.5;
        mu = 1; phi = 1; thickness = 1; p_I = 2; p_0 = 1;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
    case 2 % fork
        length_scale = 0.1; var_matern = 1; nu_matern = 2.5;
        mu = 0.25; phi = 0.75; thickness = 1; p_I = 100; p_0 = 0;
        observation_times = 1.5*linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 1.5*0.92^2*mu*phi/(2*(p_I-p_0));
    case 3 % annulus
        length_scale = 0.2; var_matern = 0.25; nu_matern = 1.5;
        mu = 0.5; phi = 0.5; thickness = 1; p_I = 6; p_0 = 1;
        observation_times = linspace(0.2,0.9,5).^2*mu*phi/(2*(p_I-p_0));
        T = 0.92^2*mu*phi/(2*(p_I-p_0));
end

matern_args = [var_matern,length_scale,nu_matern];

end