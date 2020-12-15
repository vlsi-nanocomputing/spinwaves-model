function [out_S,out_C] = DC2_ver3(in_signal)

% This function describes the behavior of the ideal DC2 (without damping).
% It receives a signal, and gives 2 output signals(out_S,out_C).
% The function has one constraint:
% *) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad]]

cd common
SW_parameters % script
cd ..

%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;        % length of the coupling region  [nm]
L_DC2=L2;       % length of the DC1 [nm], we currently use the same value of L2
gap2=10;        % the gap between the second coupled waveguides [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
limitation=0.74;

DC2_design = [h, w, d, B];
cd common
[wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC2 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC2_ff1=wm1./(2*pi);
DC2_ff2=wm2./(2*pi);


DC2_akx = in_signal(1);
dx = L2/100;
delta_phase = 0;
for i1=dx:dx:L2
    DC2_akx = DC2_akx*exp(-dx/x_freepath);
    DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(DC2_akx).^2;
    DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(DC2_ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(DC2_ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
end
DC2_Lc_avg = pi*L2/delta_phase; % average Lc
DC2_pow_par = cos(pi*L2/(2*DC2_Lc_avg))^2; % [%]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep the same frequency and phase
out_S(2) = in_signal(2);
out_C(2) = in_signal(2);
out_S(3) = in_signal(3);
out_C(3) = in_signal(3);

out_S(1) = in_signal(1) * sqrt(DC2_pow_par);
out_C(1) = in_signal(1) * sqrt(1-DC2_pow_par);


out_S(4) = in_signal(4) + tpd_DC2;
out_C(4) = in_signal(4) + tpd_DC2;


% amplitude reduction due to the damping
out_S(1) = out_S(1) * exp(-L_DC2/x_freepath);
out_C(1) = out_C(1) * exp(-L_DC2/x_freepath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

