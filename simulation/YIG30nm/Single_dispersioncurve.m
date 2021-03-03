clear
close all

h=10;                                    %thinckness (nm)
w=30;                                    %width(nm)

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
kmax=0.02;

limiation=10;

for kx=dkx:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral
       

    wm(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));
    dwm_dwh(i1)=(Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+2.*Wh)./(2.*wm(i1));
       
    i1=i1+1;

end
    

 k1=dkx:dkx:kmax;
 ff=wm./(2*pi);
 
 dw=diff(wm);
 dk=diff(k1);
 N_length=length(dk);
 
 velocity=abs(dw./dk/1000);
 lifetime=1./((damping+r*dH0*1e-9./(2*wm)).*wm.*dwm_dwh)';
 distance=abs(lifetime(1:N_length).*velocity');
 

figure
[ax1]=plotyy(k1,ff,k1(1:N_length),velocity);
xlabel('Wavenumber k (rad/nm^{-1})')
ylabel(ax1(1),'Frequency (GHz)')
ylabel(ax1(2),'Group Velocity (km/s)')

figure
[ax2]=plotyy(k1(1:N_length),lifetime(1:N_length),k1(1:N_length),distance(1:N_length));
xlabel('Wavenumber k (rad/nm^{-1})')
ylabel(ax2(1),'Lifetime (ns)')
ylabel(ax2(2),'Decay length (um)')

% interp1(k1(1:N_length),distance(1:N_length),interp1(ff,k1,2.39))
% group velocity: 206 m/s
