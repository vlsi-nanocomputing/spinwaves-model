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

%%%%%%%%%%%%%%%%%% physical parameters (constants) %%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;       % Ms  [A/m]
A=3.5e-12;      % exchange constant  [J/m]
u=pi*4e-7;      % permeability of vacuum  [H/m]
r=2.21e5;       % gyromagnetic ratio  [m/(s.A)]

j=1;

H=B*796;        % external field  [A/m] 
                % H=B/u0, u0=4*pi*e-7  [H/m]

damping=2e-4;   % damping
dH0=0.2*796;    % inhomogeneous linewidth
    
Wm=r*Ms*1e-9;   % [GHz]
Wh=r*H*1e-9;    % [GHz]
Le=sqrt(2*A/(u*Ms^2))*1e9;   % exchange length

weff=inf;       % effective width extract from Mumax3

k=1*pi/weff;
Ky=k;           %the effective wave number describing SW mode across the width direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
kmin=0.001;
limitation=8;
i1=1;

for kx=dkx:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limitation,limitation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limitation,limitation);                                                                                                       %integral
    Fkxyyd(i1)=integral(f3,-limitation,limitation);                                                                                                       %integral
    Fkxzzd(i1)=integral(f4,-limitation,limitation); 

    
    wm0(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));
    wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    
    
    % nonlinear shift coefficient Tkx calculation
    Akx(i1) = Wh+Wm/2.*(2*(Le.^2).*(kx.^2+Ky.^2)+Fkxyy0(i1)+Fkxzz0(i1));  % [rad/ns]
    Bkx(i1) = Wm/2.*(Fkxyy0(i1)-Fkxzz0(i1));  % [rad/ns]

    f5=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*(2*kx).^2./((2*kx).^2+ky.^2).*(1-(1-exp(-sqrt((2*kx).^2+ky.^2).*h))./(sqrt((2*kx).^2+ky.^2).*h))./(2*pi);
    F2kxxx0(i1) = integral(f5,-limitation,limitation); 
    
    DC1_Tkx(i1) = (Wh-Akx(i1)+Bkx(i1).^2./(2*(2*pi*in_A(2)).^2).*( Wm.*(4*Le.^2.*(kx.^2+Ky.^2)+F2kxxx0(i1))+3*Wh))./(2*pi);  % [GHz]
    
    
    i1=i1+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC1 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC1_ff0=wm0./(2*pi);
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
    display("DC1 unpredicted case: the input signals have not the same frequency")
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
elseif (in_B(3)-in_A(3)) ~= pi/2
    display("DC1 unpredicted case: the input B is not shifted by pi/2 with respect to the intput A")
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

