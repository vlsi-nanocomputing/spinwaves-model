% half-adder(HA)
if model == 2
    DCs_YIG100nm
else
    DCs_YIG30nm
end

% Area
area_HA = area_DC1 + area_DC2 + area_regS + area_regC; % [um^2]


% Energy
energy_HA = 0.5*(energy_sw + energy_sw) + energy_amp_regS + energy_amp_regC;      % [aJ]


% Propagation delay
tpd_HA_S = tpd_DC1 + tpd_DC2 + tpd_regS;  % [ns]
tpd_HA_C = tpd_DC1 + tpd_DC2 + tpd_regC;  % [ns]
tpd_HA = max(tpd_HA_S, tpd_HA_C);  % [ns]
