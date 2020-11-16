function [output] = NOT_ver1(input)


% It is an HA with one input fixed at the 1 


[S,C] = HA_ver1(input,1);
output=S;


end

