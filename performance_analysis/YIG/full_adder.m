% full-adder (FA)
half_adder

% Area
area_FA = 3*area_HA;

% Energy consumption
energy_FA = 0.5*(3*energy_sw) + 3*energy_amp_regS + 2*energy_amp_regC;

% Propagation delay
tpd_FA_C = 2*tpd_HA_S + tpd_HA_C;
tpd_FA_S = 2*tpd_HA_S;
tpd_FA = max(tpd_FA_C, tpd_FA_S);
