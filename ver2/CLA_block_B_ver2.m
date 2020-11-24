function [Ci_out, Cj1, Gik, Pik] = CLA_block_B_ver2(Gj1_k, Pj1_k, Gij, Pij, Ci)

% FAN-IN 2, duplication of the signals
[Pj1_k_1, Pj1_k_2] = duplicator_ver2(Pj1_k);
[Gij_1, Gij_2] = duplicator_ver2(Gij);
[Pij_1, Pij_2] = duplicator_ver2(Pij);
[Ci_1, Ci_2] = duplicator_ver2(Ci);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% outputs generation
Ci_out = Ci_1;
Cj1 = OR_ver2( AND_ver2(Ci_2, Pij_1), Gij_1);
Gil = OR_ver2( AND_ver2(Gij_2, Pj1_k_1), Gj1_k);
Pik = AND_ver2(Pj1_k_2, Pij_2);

end

