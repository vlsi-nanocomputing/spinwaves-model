function [output] = NOT(input)

% It is an HA with one input fixed at the 1 
fixed_input = 1;
fixed_input = DAC(fixed_input);

[output,C] = HA(input,fixed_input);

end

