%%%%%%%%%%%%%%%%%%% antiparallel static magnetization %%%%%%%%%%%%%%%%%%

%%
clear all
close all
clc

h=30;                                    %thinckness (nm)
w=100;                                    %width(nm)
SW_frequency=2.282; % frequency of the SW, it is constant for all device [GHz]

Ms=1.4e5;                                 %Ms(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))

j=1;

H=0*79.6;                              %external field(A/m)

damping=2e-4;                          %damping
dH0=0.2*796;                           %inhomogeneous linewidth
    
Wm=r*Ms*1e-9;                             %unit(GHz)
Wh=r*H*1e-9;                              %unit(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;                %exchange length

weff=inf;                                %effective width extract from Mumax3

k=1*pi/weff;
Ky=k;                                    %the effective wave number describing SW mode across the width direction
i1=1;

dkx=1e-5;
kmax=0.03;

delta1=50; % the gap between the first coupled waveguides [nm]
delta2=10; % the gap between the second coupled waveguides [nm]
% d=w+delta1;
d=w+delta2;

limiation=0.74;

for kx=dkx:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral
    Fkxyyd(i1)=integral(f3,-limiation,limiation);                                                                                                       %integral
    Fkxzzd(i1)=integral(f4,-limiation,limiation); 

    
    wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    
    ap_wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    ap_wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    


    i1=i1+1;

end
    

 k1=dkx:dkx:kmax;
 ff1=wm1./(2*pi);
 ff2=wm2./(2*pi);
 ap_ff1=ap_wm1./(2*pi);
 ap_ff2=ap_wm2./(2*pi);
 
 
 figure
 hold on
 plot(k1,real(ff1),'LineWidth',2)
 plot(k1,real(ff2),'LineWidth',2)
 plot(k1,real(ap_ff1),'LineWidth',2)
 plot(k1,real(ap_ff2),'LineWidth',2)
 grid on
 legend('Parallel configuration symmetric mode [GHz]','Parallel configuration antisymmetric mode [GHz]','Antiparallel configuration symmetric mode [GHz]','Antiparallel configuration antisymmetric mode [GHz]') 
 hold off
 xlabel('Wavenumber k [rad/nm]','FontSize',20)
 title('Parallel configuration VS Antiparallel configuration','FontSize',20)
 
 



 
 
 
 
 
 
 
 
 
 
 
 
 %% Lc(gap)
clear all
close all
clc

h = 50;                                    %thinckness (nm)
w = 100;                                    %width(nm)
SW_kx = 0.02872; % kx of initial SW (ran/nm) 
% gap = linspace(10,100); % the gap between the coupled waveguides [nm]
% d = w + gap;
B=10; % external field [mT]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;                                 %Ms(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))

j=1;

H=B*796;     %external field(A/m), H=B/u0, the B is the excitation field [mT]
                                    % u0=4*pi*e-7 H/m

damping=2e-4;                          %damping
dH0=0.2*796;                           %inhomogeneous linewidth
    
Wm=r*Ms*1e-9;                             %unit(GHz)
Wh=r*H*1e-9;                              %unit(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;                %exchange length

weff=inf;                                %effective width extract from Mumax3

k=1*pi/weff;
Ky=k;                                    %the effective wave number describing SW mode across the width direction



limiation=0.74;
kx = SW_kx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
Fkxyy0=integral(f1,-limiation,limiation);                                                                                                       %integral
Fkxzz0=integral(f2,-limiation,limiation);  
wm = sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0)));
ff0 = wm/(2*pi);



                                                                                          
Fkxyyd = zeros(1,181);                                                                                                   
Fkxzzd = zeros(1,181);
wm1 = zeros(1,181);
wm2 = zeros(1,181);
ap_wm1 = zeros(1,181);
ap_wm2 = zeros(1,181);

i1=1;
for gap=10:0.5:100                       

    d = w + gap;
    f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyyd(i1)=integral(f3,-limiation,limiation);                                                                                                       %integral
    Fkxzzd(i1)=integral(f4,-limiation,limiation); 

    
    wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0)+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0)+Wm*Fkxzzd(i1)));
    wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0)-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0)-Wm*Fkxzzd(i1)));
    
    ap_wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0)+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0)-Wm*Fkxzzd(i1)));
    ap_wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0)-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0)+Wm*Fkxzzd(i1)));
 

    i1=i1+1;

end
    

gap = 10:0.5:100;  % (nm)
group_velocity = 353.5; % (m/s), from single dispersion
ff1 = wm1./(2*pi); % Symmetric mode (GHz)
ff2 = wm2./(2*pi); % Antisymmetric mode  (GHz)
ap_ff1 = ap_wm1./(2*pi); 
ap_ff2 = ap_wm2./(2*pi); 

delta_ff = abs(ff1-ff2);
ap_delta_ff = abs(ap_ff1-ap_ff2);
Lc = (group_velocity./(2*delta_ff))*0.001;  % coupling length (um)
ap_Lc = (group_velocity./(2*ap_delta_ff))*0.001;


figure
semilogy(gap,Lc,gap,ap_Lc,'LineWidth',2);
xlabel('Gap  [nm]','FontSize',20)
legend('Parallel configuration','Antiparallel configuration')
grid on
title('Coupling length  L_c  [um]','FontSize',20)















%%  Lc(SW_kx)
clear all
close all
clc

h = 50;                                    %thinckness (nm)
w = 100;                                    %width(nm)
gap = 10; % the gap between the coupled waveguides [nm]
d = w + gap;
B = 10; % external field [mT]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;                                 %Ms(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))

j=1;

H=B*796;     %external field(A/m), H=B/u0, the B is the excitation field [mT]
                                    % u0=4*pi*e-7 H/m

damping=2e-4;                          %damping
dH0=0.2*796;                           %inhomogeneous linewidth
    
Wm=r*Ms*1e-9;                             %unit(GHz)
Wh=r*H*1e-9;                              %unit(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;                %exchange length

weff=inf;                                %effective width extract from Mumax3

k=1*pi/weff;
Ky=k;                                    %the effective wave number describing SW mode across the width direction



dkx=1e-5;
kmax=0.055;
kmin=0.02;
limiation=0.74;

k1 = kmin:dkx:kmax;
points = size(k1);
points = points(2);

Fkxyy0 = zeros(1,points);                                                                                                       %integral
Fkxzz0 = zeros(1,points);
Fkxyyd = zeros(1,points);                                                                                                   
Fkxzzd = zeros(1,points);
wm1 = zeros(1,points);
wm2 = zeros(1,points);
ap_wm1 = zeros(1,points);
ap_wm2 = zeros(1,points);


i1=1;
for kx=kmin:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral
    Fkxyyd(i1)=integral(f3,-limiation,limiation);                                                                                                       %integral
    Fkxzzd(i1)=integral(f4,-limiation,limiation); 

    wm0(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));
    wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    
    ap_wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    ap_wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    


    i1=i1+1;

end
      

ff0=wm0./(2*pi);
ff1 = wm1./(2*pi);
ff2 = wm2./(2*pi);
delta_kx = zeros(1,points);
ap_ff1 = ap_wm1./(2*pi);
ap_ff2 = ap_wm2./(2*pi);
ap_delta_kx = zeros(1,points);

for i=1:points
    SW_frequency = interp1(k1,ff0,k1(i));
    SW_ks = interp1(abs(ff1),k1,SW_frequency);
    SW_kas = interp1(abs(ff2),k1,SW_frequency);
    delta_kx(i) = abs(SW_ks-SW_kas);
    
    ap_SW_ks = interp1(abs(ap_ff1),k1,SW_frequency);
    ap_SW_kas = interp1(abs(ap_ff2),k1,SW_frequency);
    ap_delta_kx(i) = abs(ap_SW_ks-ap_SW_kas);
end

Lc = 0.001*pi./delta_kx;
ap_Lc = 0.001*pi./ap_delta_kx;

figure
hold on
plot(k1,Lc,'LineWidth',2)
plot(k1,ap_Lc,'LineWidth',2)
hold off
grid on
xlabel('Wavenumber k_x  [rad/nm]','FontSize',20)
legend('Parallel configuration','Antiparallel configuration')
title('Coupling length  L_c  [um]','FontSize',20)
axis([0.03 0.055 0 3])


















%% Lc(thickness h)
clear all
close all
clc

% h = 50;                                    %thinckness (nm)
w = 100;                                    %width(nm)
SW_kx = 0.02872; % kx of initial SW (ran/nm) 
gap = 10; % the gap between the coupled waveguides [nm]
d = w + gap;
B = 10; % external field [mT]




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;                                 %Ms(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))

j=1;

H=B*796;     %external field(A/m), H=B/u0, the B is the excitation field [mT]
                                    % u0=4*pi*e-7 H/m

damping=2e-4;                          %damping
dH0=0.2*796;                           %inhomogeneous linewidth
    
Wm=r*Ms*1e-9;                             %unit(GHz)
Wh=r*H*1e-9;                              %unit(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;                %exchange length

weff=inf;                                %effective width extract from Mumax3

k=1*pi/weff;
Ky=k;                                    %the effective wave number describing SW mode across the width direction



limiation=0.74;
kx = SW_kx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmax = 50;  %nm
hmin = 10;  %nm
points = 1000;
step = (hmax-hmin)/(points-1);


Fkxyy0 = zeros(1,points);                                                                                                       %integral
Fkxzz0 = zeros(1,points);
Fkxyyd = zeros(1,points);                                                                                                   
Fkxzzd = zeros(1,points);
wm1 = zeros(1,points);
wm2 = zeros(1,points);
ap_wm1 = zeros(1,points);
ap_wm2 = zeros(1,points);
wm0 = zeros(1,points);
Vgr = zeros(1,points);


i1=1;

for h=hmin:step:hmax
    
    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f3=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f4=@(ky)exp(i*ky.*d).*(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral
    Fkxyyd(i1)=integral(f3,-limiation,limiation);                                                                                                       %integral
    Fkxzzd(i1)=integral(f4,-limiation,limiation); 

    
    wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));
    wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    ap_wm1(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))-Wm*Fkxzzd(i1)));
    ap_wm2(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))-Wm*Fkxyyd(i1)).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+Wm*Fkxzzd(i1)));


    i1=i1+1;
end




i1=1;
for h=hmin:step:hmax
    % to calculate the group velocity, we need some values around the
    % SW_kx
    i=1;
    for kx=(SW_kx-0.001):0.0001:(SW_kx+0.001)

        f1_kx=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
        f2_kx=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

        Fkxyy0_kx(i)=integral(f1_kx,-limiation,limiation);                                                                                                       %integral
        Fkxzz0_kx(i)=integral(f2_kx,-limiation,limiation);                                                                                                       %integral
        wm(i)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0_kx(i))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0_kx(i))));

        i=i+1;
    end
    k1 = [(SW_kx-0.001):0.0001:(SW_kx+0.001)];
    dw = diff(wm);
    dk = diff(k1);
    Vgr_vector = abs(dw./dk); % m/s
    Vgr(i1) = interp1(k1(1:end-1),Vgr_vector,SW_kx);  % group velocity for fixed h and SW_kx
    i1=i1+1;
end






h = linspace(hmin,hmax,points);  %nm
% ff0=wm0./(2*pi);
ff1 = wm1./(2*pi);
ff2 = wm2./(2*pi);  
ap_ff1 = ap_wm1./(2*pi);
ap_ff2 = ap_wm2./(2*pi);  

delta_ff = abs(ff1-ff2);  % GHz
ap_delta_ff = abs(ap_ff1-ap_ff2);  % GHz


Lc = 0.001*Vgr./(2*delta_ff); % um 
ap_Lc = 0.001*Vgr./(2*ap_delta_ff); % um 


figure
hold on
plot(h,Lc,'LineWidth',2)
plot(h,ap_Lc,'LineWidth',2)
hold off
grid on
xlabel('Thickness h  [nm]','FontSize',20)
legend('Parallel configuration','Antiparallel configuration')
title('Coupling length  L_c  [um]','FontSize',20)

