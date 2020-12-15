% half-adder(HA)
directional_couplers

% Area
A_HA = A_DC1 + A_DC2 + A_regS + A_regC; % [um^2]

% Energy
E_HA = E_a_b + E_amp_regS;      % [aJ]

% Propagation delay
% tpd_HA_S = tpd_DC1 + tpd_DC2 + tpd_regS;  % [ns]
% tpd_HA_C = tpd_DC1 + tpd_DC2 + tpd_regC;  % [ns]
% tpd_HA = max(tpd_HA_S, tpd_HA_C);  % [ns]

tpd_HA_S = 150;  % [ns]
tpd_HA_C = 150;  % [ns]
tpd_HA = 150;  % [ns]