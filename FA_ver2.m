function [S,C] = FA_ver2(A,B,carry_in)


[S1,C1] = HA_ver2(A,B);

[S,C2] = HA_ver2(S1,carry_in);

C = XOR_ver2(C1,C2);
end

