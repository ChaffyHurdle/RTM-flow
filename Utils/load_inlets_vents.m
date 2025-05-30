function [inlet_func,vent_func,polar] = load_inlets_vents(case_number)

% Code inlet/vent positions
switch case_number
    case 1
        inlet_func = @(x) x(1) == 0;
        vent_func = @(x) x(1) == 1;
        polar = false;

    case 2
        inlet_func = @(x) x(1) == 0 && x(2) < 0.05;
        vent_func = @(x) x(1) == 1 && x(2) > 0.95;
        polar = true;
    case 3
        inlet_func = @(x) x(1) == 0;
        vent_func = @(x) x(1) == 1;
        polar = false;
    case 4
        r1 = 0.25;
        r2 = 1.0;
        inlet_func = @(x) abs(x(1)^2 + x(2)^2 - r1^2) < 0.002;
        vent_func = @(x) abs(x(1)^2 + x(2)^2 - r2^2) < 0.002;
        polar = true;
end

end