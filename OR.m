function [OR_out] = OR(in_A,in_B)

[S,C] = HA(in_A,in_B);
[S,C] = HA(S,C);

OR_out = S;
end

