function [OR_out] = OR_ver2(in_A,in_B)

[S,C] = HA_ver2(in_A,in_B);
[S,C] = HA_ver2(S,C);

OR_out = S;
end
