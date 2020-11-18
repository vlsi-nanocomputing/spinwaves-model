function [X] = RCA_2bit_ver2(A,B,carry_in)

X = zeros(3,3);
[X(3,:),C] = FA_ver2(A(2,:),B(2,:),carry_in);

[X(2,:),X(1,:)] = FA_ver2(A(1,:),B(1,:),C);


end

