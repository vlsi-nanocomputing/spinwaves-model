function [output] = AND4_ver2(A,B,C,D)
% AND port with 4 inputs

AND_out1 = AND_ver2(A,B);
AND_out2 = AND_ver2(C,D);

output = AND_ver2(AND_out1, AND_out2);
end

