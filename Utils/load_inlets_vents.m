function [inlet_func,vent_func] = load_inlets_vents(case_number)

% Inlet/vent functions for each case
switch case_number
    case 1 % square
        inlet_func = @(x) x(1) == 0;
        vent_func = @(x) x(1) == 1;
    case 2 % fork
        inlet_func = @(x) x(1) == 0;
        vent_func = @(x) x(1) == 1;
    case 3 % annulus
        r1 = 0.25;
        r2 = 1.0;
        inlet_func = @(x) abs(x(1)^2 + x(2)^2 - r1^2) < 0.002;
        vent_func = @(x) abs(x(1)^2 + x(2)^2 - r2^2) < 0.002;
end

end