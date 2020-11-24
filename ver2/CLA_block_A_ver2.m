function [gi,pi,si] = CLA_block_A_ver2(ai,bi,ci)
% block A of the Carry-Lookahead Adder

[HA1_S,HA1_C] = HA_ver2(ai,bi);

[gi, HA2_A] = duplicator_ver2(HA1_C);
[HA2_B, HA3_A]= duplicator_ver2(HA1_S);

[pi, xx] = HA_ver2(HA2_A,HA2_B); 

[si, xx] = HA_ver2(HA3_A,ci);

end

