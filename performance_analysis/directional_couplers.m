

%%%%%%%%%%%%%%%%%%%%%%%%% Spin-wave parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
SW_frequency = 2.282;  % GHz
SW_amplitude = 0.153;  % dimensionless

N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% design parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
vgr = 195.3;    % group velocity [m/s]
L1 = 370;       % length of coupling region  [nm]
L2 = 3000;      % length of coupling region  [nm]
L_S = 986;      % length of coupling region, regS = regenerator_S  [nm]
L_C = 566;      % length of coupling region, regC = regenerator_C  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% DCs length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_DC1 = L1+2*(5*h)/sin(0.3491);     % length of DC1  [nm]
L_DC2 = L2+2*(5*h)/sin(0.3491);     % length of DC2  [nm]
L_regS = L_S+2*(5*h)/sin(0.3491);  % length of regenerator_S  [nm]
L_regC = L_C+2*(5*h)/sin(0.3491);  % length of regenerator_C  [nm]
L_dup = L1+2*(5*h)/sin(0.3491);     % length of duplicator  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% DCs width %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_DC1 = 2*w+4*5*h;      % width of DC1  [nm]
w_DC2 = 2*w+4*5*h;      % width of DC2  [nm]
w_regS = 2*w+4*5*h;     % width of regenerator_S  [nm]
w_regC = 2*w+4*5*h;     % width of regenerator_C  [nm]
w_dup = 2*w+4*5*h;      % width of duplicator  [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Propagation delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tpd_DC1 = L_DC1/vgr;        % [ns] 
tpd_DC2 = L_DC2/vgr;        % [ns]
tpd_regS = L_regS/vgr;      % [ns]
tpd_regC = L_regC/vgr;      % [ns]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_DC1 = L_DC1/1000 * w_DC1/1000;      % [um^2]
A_DC2 = L_DC2/1000 * w_DC2/1000;      % [um^2]
A_regS = L_regS/1000 * w_regS/1000;   % [um^2]
A_regC = L_regC/1000 * w_regC/1000;   % [um^2]
A_dup = L_dup/1000 * w_dup/1000;      % [um^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy consumption %%%%%%%%%%%%%%%%%%%%%%%%%%%
E_amp_regS = 3;   % energy consumption of the amplifier of regS [aJ]
E_amp_dup = 1.5;  % energy consumption of the amplifier of duplicatore [aJ]
E_a_b = 24.6;     % energy consumption to excite inputs of HA (a,b) [aJ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


