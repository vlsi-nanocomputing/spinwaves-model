function [output] = NOT_ver2(input)

% It is an HA with one input fixed at the 1 
fixed_input = [0.153, 2.282, 0];

[S,C] = HA_ver2(input,fixed_input);
output=S;



end

