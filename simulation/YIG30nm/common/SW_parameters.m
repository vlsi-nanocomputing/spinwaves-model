% SW_parameters
% this script contains all information about the spin-wave

SW_frequency = 2.29;  % GHz
SW_amplitude = 0.093;  % dimensionless


x_freepath = 9.07e3;   % decay length [nm], calculated from 
                       % Single_dispersioncurve script 
limitation1 = 10;
limitation2 = 10;
dx = 50; % nm

N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]


%%%%%%%%%%%%%%%%%%%%%%%%%%% delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;           % thinckness  [nm]
vgr = 205.7;    % group velocity [m/s]
L1=205;         % length of the coupling region  [nm]
L2=2460;        % length of the coupling region  [nm]
L_S = 986;      % regenerator_S [nm]
L_C = 577;      % regenerator_S [nm]



L_DC1 = L1+2*(5*h)/sin(0.3491);     % length of DC1 [nm]
L_DC2 = L2+2*(5*h)/sin(0.3491);     % length of DC2 [nm]
L_regS =  L_S+2*(5*h)/sin(0.3491);  % length of regenerator_S [nm]
L_regC =  L_C+2*(5*h)/sin(0.3491);  % length of regenerator_C [nm]

% propagation delay
tpd_DC1 = L_DC1/vgr;        % [ns] 
tpd_DC2 = L_DC2/vgr;        % [ns]
tpd_regS = L_regS/vgr;      % [ns]
tpd_regC = L_regC/vgr;      % [ns]



