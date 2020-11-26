function [gi, p_i, ai_xor_bi] = CLA_block_A1_ver2(ai, bi)
% first part of block A (forth) of the Carry-Lookahead Adder 

[HA1_S,HA1_C] = HA_ver2(ai,bi);

[gi, HA2_A] = duplicator_ver2(HA1_C);
[HA2_B, ai_xor_bi]= duplicator_ver2(HA1_S);

[p_i, xx] = HA_ver2(HA2_A,HA2_B); 


end

