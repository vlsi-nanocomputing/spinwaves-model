function [output] = inverter(input)


% It is an HA with one input fixed at the 1 


[S,C] = HA(input,1);
output=S;


end

