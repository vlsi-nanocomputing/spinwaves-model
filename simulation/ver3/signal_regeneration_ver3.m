clear all
close all
clc
% This script is used to find the length (Lw) of the regenerator design.
% At the beginning, you need to set your degrated values. The goal is to
% attenuate (using the nonlinearity of DC) the degrated '0' towards zero,
% but keep the level of the '1'. At the end we can use a small amplifier
% to amplify the '1'.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S_10 = 0.032112169168809*3; % output 'S'
% S_01 = 0.032656480491817*3;
% S_11 = 0.010824493934444*3;

% dopo 1548nm
% S_10 = 0.075895597436326; % output 'S'
% S_01 = 0.078135590777412;
% S_11 = 8.970318279303576e-05;

% dopo 1107nm
S_10 = 0.066943203430362*sqrt(2); % output 'S'
S_01 = 0.069520336562717*sqrt(2);
S_11 = (4.754422686032752e-05)*sqrt(2);


% S_10 = 0.003107657797394; % C
% S_01 = 0.002429280181317;
% S_11 = 0.064074398667652;
h=30;                                    %thinckness (nm)
w=100;                                    %width(nm)
SW_frequency=2.282; % frequency of the SW, it is constant for all device [GHz]
gap=10; % the gap between the second coupled waveguides [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=w+gap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;                                 %Ms(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))

j=1;

B=0; %  field [mT]
H=B*796;     %external field(A/m), H=B/u0, u0=4*pi*e-7 H/m

damping=2e-4;                          %damping
dH0=0.2*796;                           %inhomogeneous linewidth
    
Wm=r*Ms*1e-9;                             %unit(GHz)
Wh=r*H*1e-9;                              %unit(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;                %exchange length

weff=inf;                                %effective width extract from Mumax3

k=1*pi/weff;
Ky=k;                                    %the effective wave number describing SW mode across the width direction
i1=1;

dkx=1e-3;
kmax=0.04;
kmin=0.001;



limitation=0.63;

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
    
    Tkx(i1) = (Wh-Akx(i1)+Bkx(i1).^2./(2*wm0(i1).^2).*( Wm.*(4*Le.^2.*(kx.^2+Ky.^2)+F2kxxx0(i1))+3*Wh))./(2*pi);  % [GHz]
    
    i1=i1+1;
end
    

k1=dkx:dkx:kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);



delta_ph10=0;
delta_ph01=0;
delta_ph11=0;
x_freepath = 8.58e3;
dl=1;
Lw=1:1:903;

for i1=1:1:903
    S_10 = S_10*exp(-dl/x_freepath);
    S_01 = S_01*exp(-dl/x_freepath);
    S_11 = S_11*exp(-dl/x_freepath);

    ff1_s10 = ff1+Tkx.*abs(S_10).^2;
    ff2_s10 = ff2+Tkx.*abs(S_10).^2;
    
    ff1_s01 = ff1+Tkx.*abs(S_01).^2;
    ff2_s01 = ff2+Tkx.*abs(S_01).^2;
    
    ff1_s11 = ff1+Tkx.*abs(S_11).^2;
    ff2_s11 = ff2+Tkx.*abs(S_11).^2;
    
    ks = interp1(abs(ff1_s10),k1,SW_frequency);
    kas = interp1(abs(ff2_s10),k1,SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph10 = delta_ph10 + delta_k*dl;
    
    ks = interp1(abs(ff1_s01),k1,SW_frequency);
    kas = interp1(abs(ff2_s01),k1,SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph01 = delta_ph01 + delta_k*dl;
    
    ks = interp1(abs(ff1_s11),k1,SW_frequency);
    kas = interp1(abs(ff2_s11),k1,SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph11 = delta_ph11 + delta_k*dl;

    
    pow_par10(i1) = cos(delta_ph10/2)^2;
    pow_par01(i1) = cos(delta_ph01/2)^2;
    pow_par11(i1) = cos(delta_ph11/2)^2;
    
end


out10 = S_10*sqrt((pow_par10(end)));
out01 = S_01*sqrt((pow_par01(end)));
out11 = S_11*sqrt((pow_par11(end)));
% normalization([out10,out01,out11])

x = linspace(0,1,1000);
y = x*(S_10/S_01)^2;
figure
plot(x,y)

figure
hold on
plot(Lw,pow_par10,'LineWidth',1.5)
plot(Lw,pow_par01,'LineWidth',1.5)
plot(Lw,pow_par11,'LineWidth',1.5)
hold off
xlabel('L_w_3  [nm]','FontSize',20)
legend('10','01','11')
% clegend('10','01')
title('Normalized output power','FontSize',15)