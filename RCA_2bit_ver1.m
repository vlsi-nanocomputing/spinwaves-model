function [X] = RCA_2bit_ver1(A,B,carry_in)

X = zeros(1,3);

[X(3),C] = FA_ver1(A(2),B(2),carry_in);
[X(2),X(1)] = FA_ver1(A(1),B(1),C);


end

