% full-adder (FA)
half_adder

% Area
A_FA = 3*A_HA;

% Energy consumption
E_FA = 3*E_a_b + 3*E_amp_regS;

% Propagation delay
tpd_FA_C = 2*tpd_HA_S + tpd_HA_C;
tpd_FA_S = 2*tpd_HA_S;
tpd_FA = max(tpd_FA_C, tpd_FA_S);
