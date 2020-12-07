function [mux2to1_out] = mux2to1_ver2(in_0,in_1,sel)

% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 

[sel1,sel2] = duplicator_ver2(sel);

out1 = AND_ver2(in_0,NOT_ver2(sel1));
out2 = AND_ver2(in_1,sel2);

mux2to1_out = OR_ver2(out1,out2);



end

