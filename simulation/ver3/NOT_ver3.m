function [output] = NOT_ver3(input)

% It is an HA with one input fixed at the 1 
fixed_input = 1;
fixed_input = DAC_ver3(fixed_input);

[output,C] = HA_ver3(input,fixed_input);

end

