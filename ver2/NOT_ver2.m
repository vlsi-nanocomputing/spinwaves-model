function [output] = NOT_ver2(input)

% It is an HA with one input fixed at the 1 
fixed_input = 1;
fixed_input = DAC_ver2(fixed_input);

[S,C] = HA_ver2(input,fixed_input);
output=S;



end

