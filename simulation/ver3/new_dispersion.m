clear all
close all
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;         % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
DC2_akx=2*0.072290169039954/sqrt(2);  % input SW of the DC2. Microwave field = 2mT 
SW_frequency = 2.282; % GHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;


limitation=0.715;
% limitation=0.6517; % 0.6512 with damping
% limitation=0.6512; 
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
    
    Tkx(i1) = (Wh-Akx(i1)+Bkx(i1).^2./(2*wm0(i1).^2).*( Wm.*(4*Le.^2.*(kx.^2+Ky.^2)+F2kxxx0(i1))+3*Wh))./(2*pi);  % [GHz]    
    i1=i1+1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Directional coupler 2 
%    DC2 input  ----------------------------------------------- output S
%                      /                           \ 
%                     /                        'l'  \
%                                                    \
%                                                      -------- output 
%
% Now we add the contribution of the region 'l', 
% which has a length = 585 nm. 

% % % % % % % % len=585;
% % % % % % % % % Let us consider a finite partitioning of the 'len' into 100 subintervals
% % % % % % % % dl = len/100;
% % % % % % % % dgap=200/100;
% % % % % % % % i2=1;
% % % % % % % % for l=dl:dl:len
% % % % % % % %     d = d + i2*dgap;
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     limitation=0.74;
% % % % % % % %     i1=1;
% % % % % % % %     % for each dl, we calculate ff1(wm1) and ff(wm2)
% % % % % % % %     for kx=dkx:dkx:kmax                       
% % % % % % % % 
% % % % % % % %         f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
% % % % % % % %         f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
% % % % % % % %         f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
% % % % % % % %         f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
% % % % % % % % 
% % % % % % % %         Fkxyy0(i1)=integral(f1,-limitation,limitation);                                                                                                       %integral
% % % % % % % %         Fkxzz0(i1)=integral(f2,-limitation,limitation);                                                                                                       %integral
% % % % % % % %         Fkxyyd(i1)=integral(f3,-limitation,limitation);                                                                                                       %integral
% % % % % % % %         Fkxzzd(i1)=integral(f4,-limitation,limitation); 
% % % % % % % % 
% % % % % % % %         wm1_t(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
% % % % % % % %         wm2_t(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
% % % % % % % % 
% % % % % % % % 
% % % % % % % %         i1=i1+1;
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %     wm1 = wm1 + wm1_t;
% % % % % % % %     wm2 = wm2 + wm2_t;
% % % % % % % %     
% % % % % % % %     i2=i2+1;
% % % % % % % % end
% % % % % % % % 
k1=dkx:dkx:kmax;
ff1 = wm1/(2*pi);
ff2 = wm2/(2*pi);


x_freepath = 8.58e3;
dx = L2/100;
delta_phase = 0;
for i1=dx:dx:L2
    DC2_akx = DC2_akx*exp(-dx/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
end
DC2_Lc_avg = pi*L2/delta_phase; % average Lc
DC2_pow_par = cos(pi*L2/(2*DC2_Lc_avg))^2; % [%]

% % % % % % % % % nonlinear frequency shift
% % % % % % % % ff1_s = ff1 + Tkx.*abs(DC2_akx).^2; 
% % % % % % % % ff2_s = ff2 + Tkx.*abs(DC2_akx).^2;
% % % % % % % % 
% % % % % % % % k1=dkx:dkx:kmax;
% % % % % % % % % hold on
% % % % % % % % % plot(k1,ff1)
% % % % % % % % % plot(k1,ff2)
% % % % % % % % % hold off
% % % % % % % % % xlabel('Wavenumber k [rad/nm]','FontSize',20)
% % % % % % % % % legend('ff1\_s  [GHz]','ff2\_s  [GHz]')
% % % % % % % % % 
% % % % % % % % ks = interp1(real(ff1_s),k1,SW_frequency);
% % % % % % % % kas = interp1(real(ff2_s),k1,SW_frequency);
% % % % % % % % Lc = pi/abs(ks-kas);  % [nm]
% % % % % % % % % pow_par = cos(pi*L2/(2*Lc))^2; %  Pout1/(Pout1+Pout2)