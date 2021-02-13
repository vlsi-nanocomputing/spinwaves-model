function [output] = AND4_ver3(A,B,C,D)
% AND port with 4 inputs

AND_out1 = AND_ver3(A,B);
AND_out2 = AND_ver3(C,D);

output = AND_ver3(AND_out1, AND_out2);
end

