function [out_signal] = regenerator_C_ver2(in_signal)

% This function describes the behavior of the regenerator (DC+amplifier)
% fot the carry bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad]]


%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain=1.01761;   % gain of the VCMA amplifier
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
Lw=577;         % length of the coupling region  [nm]
gap=10;        % the gap between the second coupled waveguides [nm]
d=w+gap;       % [nm]
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.04;
limitation=0.74;

DC_design = [h, w, gap, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC_design);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% regenerator operation %%%%%%%%%%%%%%%%%%%%%%%%%%%

k1=dkx:dkx:kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);


% nonlinear shift due to the input power
akx = in_signal(1);
ff1_s = ff1+Tkx.*abs(akx).^2;
ff2_s = ff2+Tkx.*abs(akx).^2;

ks = interp1(abs(ff1_s),k1,in_signal(2));
kas = interp1(abs(ff2_s),k1,in_signal(2));
Lc = pi/abs(ks-kas);  % [nm]

pow_par = cos(pi*Lw/(2*Lc))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = in_signal; % to have the same frequency and phase

% propagation delay
cd common
SW_parameters % script
cd ..
out_signal(4) = in_signal(4) + tpd_regC; 

% power splitting and amplification
out_signal(1) = in_signal(1) * sqrt(pow_par);
out_signal = amplifier_ver2(out_signal,gain);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

