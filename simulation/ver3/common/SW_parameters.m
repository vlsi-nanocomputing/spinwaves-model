% SW_parameters
% this script contains all information about the spin-wave

SW_frequency = 2.282;  % GHz
% SW_amplitude = 0.0781;  % dimensionless
SW_amplitude = 0.0779;  % dimensionless
%%%%%%%%%%%%%%%%%%%%%%%% 2mT=0.051016316539556;
%%%%%%%%%%%%%%%%%%%%%%%% 4mT=0.102233737482546;



x_freepath = 8.58e3;   % decay length [nm], calculated from 
                       % Single_dispersioncurve script 
limitation1 = 0.53;
limitation2 = 0.63;
dx = 50; % nm

N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]


%%%%%%%%%%%%%%%%%%%%%%%%%%% delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
vgr = 205.7;    % group velocity [m/s]
L1=370;         % length of the coupling region  [nm]
L2=3000;        % length of the coupling region  [nm]
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



