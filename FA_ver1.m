function [S,C] = FA_ver1(A,B,carry_in)


[S1,C1] = HA_ver1(A,B);

[S,C2] = HA_ver1(S1,carry_in);

C = OR_ver1(C1,C2);
end

