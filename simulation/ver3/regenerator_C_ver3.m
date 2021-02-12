function [out_signal] = regenerator_C_ver3(in_signal)

% This function describes the behavior of the regenerator (DC+amplifier)
% fot the carry bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad], delay]


%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_out = 3.6;   % 
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L4=1950;
gap=10;        % the gap between the second coupled waveguides [nm]
d=w+gap;       % [nm]
B=0;            % external field [mT]
DC4_akx = in_signal(1);
cd common
SW_parameters % script
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.04;
limitation2=0.63;

DC_design = [h, w, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation2, DC_design);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% regenerator operation %%%%%%%%%%%%%%%%%%%%%%%%%%%

k1=dkx:dkx:kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);


delta_phase = 0;
dx = L4/65;
for i1=dx:dx:L4
    DC4_akx = DC4_akx*exp(-dx/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC4_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC4_akx).^2;
    DC4_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC4_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC4_ks-DC4_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
end


Lc_avg = pi*L4/delta_phase;  % [nm]

pow_par = cos(pi*L4/(2*Lc_avg))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = in_signal; % to have the same frequency and phase

% propagation delay
cd common
SW_parameters % script
cd ..
out_signal(4) = in_signal(4) + tpd_regC; 

% power splitting and amplification
out_signal(1) = DC4_akx * sqrt(pow_par);
out_signal = amplifier_ver3(out_signal,gain_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

