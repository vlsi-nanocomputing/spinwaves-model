
% parameters interface
% nand2
INV_NAND2_fanout1 

area_NAND2 = area_NAND2; % um^2
tpd_NAND2 = tau_NAND2;   % ns
ps_NAND2 = Ps_NAND2;     % W
pdyn_NAND2 = Pdyn_NAND2; % W*s = J



%%%%%%% FA parameters
tpd_FA_S = 6*tpd_NAND2;                     % FA sum bit delay, [ns]
tpd_FA_C = 5*tpd_NAND2;                     % FA carry bit delay, [ns]
area_FA = 9*area_NAND2;                     % FA area, [um^2]
ps_FA = 9*ps_NAND2;                         % FA static power, [W]
pdyn_FA = 9*pdyn_NAND2;                     % FA dynamic power [W*s = J]

