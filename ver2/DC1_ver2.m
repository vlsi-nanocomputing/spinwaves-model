function [out,out_I] = DC1_ver2(in_A,in_B)

% This function describes the behavior of the ideal DC1 (without damping).
% It receives 2 signals (A and B), and gives 2 output signals(out,out_I).
% The function has some constraints:
% 1) the 2 input signals have the same frequency
% 2) phase(B) - phase(A) = pi/2
% 3) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad]]



%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L1=370;         % length of the coupling region  [nm]
gap1=50;        % the gap between the coupled waveguides  [nm]
d=w+gap1;       % [nm]
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
limitation=8;

DC1_design = [h, w, gap1, d, B];
cd common
[wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC1 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC1_ff1=wm1./(2*pi);
DC1_ff2=wm2./(2*pi);

% nonlinear shift due to the input power
DC1_akx = in_A(1);
DC1_ff1_s = DC1_ff1 + DC1_Tkx.*abs(DC1_akx).^2;
DC1_ff2_s = DC1_ff2 + DC1_Tkx.*abs(DC1_akx).^2;

DC1_ks = interp1(real(DC1_ff1),k1,in_A(2));   % [rad/nm]
DC1_kas = interp1(real(DC1_ff2),k1,in_A(2));  % [rad/nm]
DC1_Lc = pi/abs(DC1_ks-DC1_kas);         % [nm]
DC1_pow_par = cos(pi*L1/(2*DC1_Lc))^2;  %  [%]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%

% the single signals (in terms of amplitude) at the 2 outputs of the DC 
% the number 1 and 2 indicate the first and second (out I) outputs,
% respectively
in_A1_akx = in_A(1)*sqrt(DC1_pow_par);
in_A2_axk = in_A(1)*sqrt(1-DC1_pow_par);
in_B1_akx = in_B(1)*sqrt(1-DC1_pow_par);
in_B2_axk = in_B(1)*sqrt(DC1_pow_par);


if in_A(2) ~= in_B(2)
    display('DC1 unpredicted case: the input signals have not the same frequency')
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
elseif (in_B(3)-in_A(3)) ~= pi/2
    display('DC1 unpredicted case: the input B is not shifted by pi/2 with respect to the intput A')
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
else
    out(1) = in_A1_akx+in_B1_akx; % constructive interference
    out(2) = in_A(2);
    out(3) = in_A(3);
    
    % destructive interference
    if in_B2_axk>in_A2_axk
        out_I(1) = in_B2_axk-in_A2_axk; % close to 0 
        out_I(2) = in_A(2);
        out_I(3) = in_B(3);  % pi/2
    else
        out_I(1) = in_A2_axk(1)-in_B2_axk(1); % close to 0
        out_I(2) = in_A(2);
        out_I(3) = in_A(3)-pi/2; % -pi/2
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













end

