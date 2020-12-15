function [S,C] = FA_ver3(A,B,carry_in)


[S1,C1] = HA_ver3(A,B);

[S,C2] = HA_ver3(S1,carry_in);

C = XOR_ver3(C1,C2);
end

