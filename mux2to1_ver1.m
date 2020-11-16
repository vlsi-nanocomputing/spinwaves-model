function [mux2to1_out] = mux2to1_ver1(in_0,in_1,sel)

% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 

out1 = AND_ver1(in_0,NOT(sel));
out2 = AND_ver1(in_1,sel);

mux2to1_out = OR_ver1(out1,out2);



end

