function [mux2to1_out] = mux2to1(in_0,in_1,sel)

% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 

out1 = AND(in_0,NOT(sel));
out2 = AND(in_1,sel);

mux2to1_out = OR(out1,out2);



end

