function [Gik, Pik] = CLA_block_B1_ver2(Gj1_k, Pj1_k, Gij, Pij)
% first part of block B (forth) of the Carry-Lookahead Adder 


[Pj1_k_1, Pj1_k_2] = duplicator_ver2(Pj1_k);

Gik = OR_ver2( AND_ver2(Gij, Pj1_k), Gj1_k);  
Pik = AND_ver2(Pj1_k_2, Pij);

end

