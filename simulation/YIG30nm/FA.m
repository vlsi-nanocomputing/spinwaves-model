function [S,C] = FA(A,B,carry_in)


[S1,C1] = HA(A,B);

[S,C2] = HA(S1,carry_in);

C = XOR(C1,C2);
end

