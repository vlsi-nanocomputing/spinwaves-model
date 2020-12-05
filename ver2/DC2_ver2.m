function [out_S,out_C] = DC2_ver2(in_signal)

% This function describes the behavior of the ideal DC2 (without damping).
% It receives a signal, and gives 2 output signals(out_S,out_C).
% The function has one constraint:
% *) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad]]


%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;        % length of the coupling region  [nm]
gap2=10;        % the gap between the second coupled waveguides [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
limitation=0.74;

DC2_design = [h, w, gap2, d, B];
cd common
[wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC2 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC2_ff1=wm1./(2*pi);
DC2_ff2=wm2./(2*pi);

% nonlinear shift due to the input power
DC2_akx = in_signal(1);
DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(DC2_akx).^2;
DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(DC2_akx).^2;

DC2_ks = interp1(real(DC2_ff1_s),k1,in_signal(2));
DC2_kas = interp1(real(DC2_ff2_s),k1,in_signal(2));
DC2_Lc = pi/abs(DC2_ks-DC2_kas);  % [nm]
DC2_pow_par = cos(pi*L2/(2*DC2_Lc))^2; % [%]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep the same frequency and phase
out_S(2) = in_signal(2);
out_C(2) = in_signal(2);
out_S(3) = in_signal(3);
out_C(3) = in_signal(3);

out_S(1) = in_signal(1) * sqrt(DC2_pow_par);
out_C(1) = in_signal(1) * sqrt(1-DC2_pow_par);

cd common
SW_parameters % script
cd ..
out_S(4) = in_signal(4) + tpd_DC2;
out_C(4) = in_signal(4) + tpd_DC2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

