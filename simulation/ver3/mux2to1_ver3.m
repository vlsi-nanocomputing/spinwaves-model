function [mux2to1_out] = mux2to1_ver3(in_0,in_1,sel)

% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 

[sel1,sel2] = duplicator_ver3(sel);

out1 = AND_ver3(in_0,NOT_ver3(sel1));
out2 = AND_ver3(in_1,sel2);

mux2to1_out = OR_ver3(out1,out2);



end

