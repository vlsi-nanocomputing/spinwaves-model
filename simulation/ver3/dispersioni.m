% clear all
% close all
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;                                    %thinckness (nm)
w=100;                                    %width(nm)
SW_frequency=2.282; % frequency of the SW, it is constant for all device [GHz]
L2=3000; %nm
gap2=10; % the gap between the second coupled waveguides [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% akx = 0.013831560061652

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=w+gap2;

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

dkx=1e-4;
kmax=0.025;
kmin=0.001;



limitation=0.74;

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
%     Akx(i1) = Wh+Wm/2.*(2*(Le.^2).*(kx.^2+Ky.^2)+Fkxyy0(i1)+Fkxzz0(i1));  % [rad/ns]
%     Bkx(i1) = Wm/2.*(Fkxyy0(i1)-Fkxzz0(i1));  % [rad/ns]
% 
%     f5=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*(2*kx).^2./((2*kx).^2+ky.^2).*(1-(1-exp(-sqrt((2*kx).^2+ky.^2).*h))./(sqrt((2*kx).^2+ky.^2).*h))./(2*pi);
%     F2kxxx0(i1) = integral(f5,-limitation,limitation); 
%     
%     DC2_Tkx(i1) = (Wh-Akx(i1)+Bkx(i1).^2./(2*wm0(i1).^2).*( Wm.*(4*Le.^2.*(kx.^2+Ky.^2)+F2kxxx0(i1))+3*Wh))./(2*pi);  % [GHz]
%     
    i1=i1+1;
end
    

k1=dkx:dkx:kmax;
DC2_ff0=wm0./(2*pi);
DC2_ff1=wm1./(2*pi);
DC2_ff2=wm2./(2*pi);

% 
% DC2_akx = 0.08053;
% akx_vec = 2*0.0138e-3/sqrt(2);
% DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(akx_vec).^2;
% DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(akx_vec).^2;
%  
% % %  
a1 = 0.0138e-3/sqrt(2);
a2 = 2*0.0138e-3/sqrt(2);
ff1_s_a1 = DC2_ff1+DC2_Tkx.*abs(a1).^2;
ff2_s_a1 = DC2_ff2+DC2_Tkx.*abs(a1).^2;
ff1_s_a2 = DC2_ff1+DC2_Tkx.*abs(a2).^2;
ff2_s_a2 = DC2_ff2+DC2_Tkx.*abs(a2).^2;
figure
hold on
plot(k1,ff1_s_a1)
plot(k1,ff2_s_a1)
plot(k1,ff1_s_a2)  % si sovrappone con ff1_s_a1 in quanto l'errore è dell'ordine del e-9
plot(k1,ff2_s_a2)
hold off
legend('ff1a1','ff2a1','ff1a2','ff2a2')



points=size(DC2_ff1);
 points=points(2);
 figure
 hold on
%  plot(k1,DC2_ff0)
 plot(k1,DC2_ff1,'LineWidth',2)
 plot(k1,DC2_ff2,'LineWidth',2)
%  plot(k1,DC2_ff1_s)
%  plot(k1,DC2_ff2_s)
 plot(k1,SW_frequency*ones(1,points),'LineWidth',2)
 grid on
% %  legend('DC2 symmetric mode [GHz]','Antisymmetric mode [GHz]','Shifted symmetric mode [GHz]','Shifted antisymmetric mode [GHz]','SW frequency (2.282 GHz)')
 legend('Symmetric mode [GHz]','Antisymmetric mode [GHz]','SW frequency (2.282 GHz)')
 hold off
 xlabel('Wavenumber k [rad/nm]','FontSize',20)
%  axis([0 0.025 1.8 2.4])

 
DC2_ks = interp1(abs(DC2_ff1_s),k1,SW_frequency);
DC2_kas = interp1(abs(DC2_ff2_s),k1,SW_frequency);
DC2_Lc = pi/abs(DC2_ks-DC2_kas);  % [nm]
DC2_pow_par = cos(pi*L2/(2*DC2_Lc))^2;

