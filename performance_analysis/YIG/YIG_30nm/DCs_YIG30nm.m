% YIG 30 nm technological parameters
% It is the lowest level

%%%%%%%%%%%%%%%%%%%%%%%%% design parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;           % thinckness  [nm]
w=30;          % width  [nm]
vgr_sing = 205.7;    % group velocity of single waveguide [m/s]
vgr_coup = 137;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
L_sing = 2*(5*h)/sin(0.3491);
L1 = 260;       % length of coupled region  [nm]
L2 = 2460;      % length of coupled region  [nm]
L_S1 = 900;      % length of couplied region, S = regenerator_S  [nm]
L_S2 = 978;      % length of coupled region, S = regenerator_S  [nm]
L_C = 1450;      % length of coupled region, C = regenerator_C  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% DCs length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_DC1 = L1 + L_sing;     % length of DC1  [nm]
L_DC2 = L2 + L_sing;     % length of DC2  [nm]
L_regS = L_S1 + L_S2;  % length of regenerator_S  [nm]
L_regC = L_C;  % length of regenerator_C  [nm]
L_dup = L_DC1;     % length of duplicator  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% DCs width %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_DC1 = 2*w+4*5*h;      % width of DC1  [nm]
w_DC2 = 2*w+4*5*h;      % width of DC2  [nm]
w_regS = 2*w+4*5*h;     % width of regenerator_S  [nm]
w_regC = 2*w+4*5*h;     % width of regenerator_C  [nm]
w_dup = 2*w+4*5*h;      % width of duplicator  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Propagation delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tpd_DC1 = L1/vgr_coup + L_sing/vgr_sing;        % [ns] 
tpd_DC2 = L2/vgr_coup + L_sing/vgr_sing;        % [ns]
tpd_regS = L_regS/vgr_coup;      % [ns]
tpd_regC = L_regC/vgr_coup;      % [ns]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_DC1 = L_DC1/1000 * w_DC1/1000;      % [um^2]
area_DC2 = L_DC2/1000 * w_DC2/1000;      % [um^2]
area_regS = L_regS/1000 * w_regS/1000;   % [um^2]
area_regC = L_regC/1000 * w_regC/1000;   % [um^2]
area_dup = L_dup/1000 * w_dup/1000;      % [um^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy consumption %%%%%%%%%%%%%%%%%%%%%%%%%%%
energy_amp_regS = 3*(10/4);   % energy consumption of the amplifier of regS [aJ]
energy_amp_regC = 3*(1.7/4);
energy_amp_dup = 3/2;  % energy consumption of the amplifier of duplicatore [aJ]
energy_sw = 1.96;     % energy consumption to excite a single input wave [aJ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






