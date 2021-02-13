function [OR_out] = OR_ver3(in_A,in_B)

[S,C] = HA_ver3(in_A,in_B);
[S,C] = HA_ver3(S,C);

OR_out = S;
end
