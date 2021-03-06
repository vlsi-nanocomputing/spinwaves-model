function [output] = AND4(A,B,C,D)
% AND port with 4 inputs

AND_out1 = AND(A,B);
AND_out2 = AND(C,D);

output = AND(AND_out1, AND_out2);
end

