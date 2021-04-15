% SW_parameters
% this script contains all information about the spin-wave

switch model 
    case 1  % YIG100nm behavioral model
        SW_frequency = 2.282;
        SW_amplitude = 0.153;
        limitation1 = 8;
        limitation2 = 0.74;
        duplicator_gain1 = 1/0.5002; 
        duplicator_gain2 = 1/(1-0.5002);
        gain_S = 2; % the amplifier gain at the DC2 out_S (without regS) 
        gain_C = 1; % the amplifier gain at the DC2 out_C (without regC)
        N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
    case 2 % YIG100nm physical model
        SW_frequency = 2.282;
        SW_amplitude = 0.0779;
        x_freepath = 8.58e3;   % decay length [nm], calculated from 
                               % Single_dispersioncurve script 
        limitation1 = 0.53;    % for gap=50nm
        limitation2 = 0.63;    % for gap=10nm
        duplicator_gain1 = 2.37;
        duplicator_gain2 = 2.3;
        gain_S = 5.5; % the amplifier gain at the DC2 out_S (without regS) 
        gain_C = 1.5; % the amplifier gain at the DC2 out_C (without regC)
        N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
    case 3 % YIG30nm physical model
        SW_frequency = 2.29;  % GHz
        SW_amplitude = 0.093;  % dimensionless
        x_freepath = 9.07e3;   % decay length [nm], calculated from 
                               % Single_dispersioncurve script 
        limitation1 = 10;
        limitation2 = 10;
        duplicator_gain1 = 2.25;
        duplicator_gain2 = 2.07;
        gain_S = 4.2; % the amplifier gain at the DC2 out_S (without regS) 
        gain_C = 1; % the amplifier gain at the DC2 out_C (without regC)
        N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
    otherwise % Logical model
        duplicator_gain1 = 2; 
        duplicator_gain2 = 2;
        gain_S = 4; % the amplifier gain at the DC2 out_S (without regS) 
        gain_C = 1; % the amplifier gain at the DC2 out_C (without regC)
        N_inf = 1;  % number of information to define a spin-wave vector
end
       

dx = 50; % nm, discretization resolution to discretize DCs
dkx=1e-3;  % to define vector kx = [dkx:dkx:kmax]
kmax=0.025;

