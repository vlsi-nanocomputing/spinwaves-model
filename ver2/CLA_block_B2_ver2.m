function [ci_out, cj1] = CLA_block_B2_ver2(ci, Gij, Pij)
% first part of block B (back) of the Carry-Lookahead Adder 

[ci1, ci2] = duplicator_ver2(ci);

ci_out = ci1;
cj1 = OR_ver2( AND_ver2(ci2,Pij), Gij);

end

