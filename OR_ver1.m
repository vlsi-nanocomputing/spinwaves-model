function [OR_out] = OR_ver1(in_A,in_B)

[S,C] = HA_ver1(in_A,in_B);
[S,C] = HA_ver1(S,C);

OR_out = S;
end

