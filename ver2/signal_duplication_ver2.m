clear all
close all
clc
% This script is used to design a block that can duplicate a signal into
% two same ones.



%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_ak=0.153;       % max amplitude
min_ak=0;           % min amplitude
h=30;               % thinckness  [nm]
w=100;              % width  [nm]
gap=50;             % the gap between the coupled waveguides  [nm]
SW_frequency=2.282; % frequency of the SW, it is constant for all device [GHz]
d=w+gap;            % [nm]
B=0;                % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
kmax=0.04;
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
    
    Tkx(i1) = (Wh-Akx(i1)+Bkx(i1).^2./(2*(2*pi*SW_frequency).^2).*( Wm.*(4*Le.^2.*(kx.^2+Ky.^2)+F2kxxx0(i1))+3*Wh))./(2*pi);  % [GHz]
    
    
    i1=i1+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1=dkx:dkx:kmax;
ff0=wm0./(2*pi);
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);


akx1 = min_ak;  % SW amplitude of the degraded '0'
ff1_s1 = ff1+Tkx.*abs(akx1).^2;
ff2_s1 = ff2+Tkx.*abs(akx1).^2;
akx2 = max_ak;  % SW amplitude of the degraded '1'
ff1_s2 = ff1+Tkx.*abs(akx2).^2;
ff2_s2 = ff2+Tkx.*abs(akx2).^2;


ks1 = interp1(abs(ff1_s1),k1,SW_frequency);
kas1 = interp1(abs(ff2_s1),k1,SW_frequency);
Lc1 = pi/abs(ks1-kas1);  % [nm]

ks2 = interp1(abs(ff1_s2),k1,SW_frequency);
kas2 = interp1(abs(ff2_s2),k1,SW_frequency);
Lc2 = pi/abs(ks2-kas2);  % [nm]

i1=1;
for Lw=100:1:3000

    pow_par1(i1) = cos(pi*Lw/(2*Lc1))^2;
    
    pow_par2(i1) = cos(pi*Lw/(2*Lc2))^2;
    
    i1=i1+1;
end

Lw=100:1:3000;
hold on
plot(Lw,pow_par1,'LineWidth',1.5)
plot(Lw,pow_par2,'LineWidth',1.5)
hold off
xlabel("Lw  [nm]")
legend('logic 0','logic 1')
