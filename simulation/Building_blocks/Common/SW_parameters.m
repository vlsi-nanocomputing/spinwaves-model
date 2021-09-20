% model_parameters
% this script contains all information about the spin-wave

switch model 
    case 'YIG 100nm'  % YIG100nm behavioral model
        model_parameters.SW_frequency = 2.282;
        model_parameters.SW_amplitude = 0.0779;
        model_parameters.vgr_sing = 195.3;    % group velocity of single waveguide [m/s]
        model_parameters.vgr_coup = 25;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
        model_parameters.x_freepath = 8.58e3;   % decay length [nm], calculated from 
                               % Single_dispersioncurve script 
        model_parameters.limitation1 = 0.53;    % for gap=50nm
        model_parameters.limitation2 = 0.63;    % for gap=10nm
        model_parameters.limitation_sing_waveg = 10; % for single waveguide
        model_parameters.duplicator_gain1 = 2.37;
        model_parameters.duplicator_gain2 = 2.3;
        model_parameters.gain_S = 5.5; % the amplifier gain at the DC2 out_S (without regS) 
        model_parameters.gain_C = 1.5; % the amplifier gain at the DC2 out_C (without regC)
        model_parameters.N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
    case 'YIG 30nm' % YIG30nm physical model
        model_parameters.SW_frequency = 2.29;  % GHz
        model_parameters.SW_amplitude = 0.093;  % dimensionless
        model_parameters.vgr_sing = 205.7;    % group velocity of single waveguide [m/s]
        model_parameters.vgr_coup = 137;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
        model_parameters.x_freepath = 9.07e3;   % decay length [nm], calculated from 
                               % Single_dispersioncurve script 
        model_parameters.limitation1 = 10;
        model_parameters.limitation2 = 10;
        model_parameters.limitation_sing_waveg = 10; % for single waveguide
        model_parameters.duplicator_gain1 = 2.25;
        model_parameters.duplicator_gain2 = 2.07;
        model_parameters.gain_S = 4.2; % the amplifier gain at the DC2 out_S (without regS) 
        model_parameters.gain_C = 1; % the amplifier gain at the DC2 out_C (without regC)
        model_parameters.N_inf = 4;  % number of information to define a spin-wave vector
            % example: [amplitude, frequency, phase, delay]
    otherwise % Logical model
        fprintf('Error, unvalid model selected')
end
       

model_parameters.dx = 50; % nm, discretization resolution to discretize DCs
model_parameters.dkx=1e-3;  % to define vector kx = [dkx:dkx:kmax]
model_parameters.kmax=0.025;

