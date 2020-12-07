function [out_S,out_C] = DC2_ver1(in_power)

% DC2 operates in a nonlinear regime and simultaneously performs XOR and AND logic operations.
% in_power = [0, 0.5, 1]
switch in_power
    case 1 
        out_S = 0;
        out_C = in_power;
    case 0.5
        out_S = in_power;
        out_C = 0;
    case 0
        out_S = 0;
        out_C = 0;
    otherwise
        disp('unpredicted case, check your DC2 input!')
        out_S = -10;
        out_C = -10;
end


end

