function [mux2to1_out] = mux2to1(in_0,in_1,sel)

% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 

[sel1,sel2] = duplicator(sel);

out1 = AND(in_0,NOT(sel1));
out2 = AND(in_1,sel2);

mux2to1_out = OR(out1,out2);



end

