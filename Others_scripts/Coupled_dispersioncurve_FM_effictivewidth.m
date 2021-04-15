clear
close all

h=30;                                     %thinckness (nm)
w=100;                                    %width(nm)
gap=10;                                   %the gap between the two waveguides (nm)
Ms=1.4e5;                                 %magnetization(A/m)
A=3.5e-12;                                %exchange constant(J/m)
u=pi*4e-7;                                %permeability of vacuum(H/m)
r=2.21e5;                                 %gyromagnetic ratio(m/(s.A))
H=0*79.6;                                 %external field(A/m)

Wm=r*Ms*1e-9;                             %units(GHz)
Wh=r*H*1e-9;                              %units(GHz)
Le=sqrt(2*A/(u*Ms^2))*1e9;
d=w+gap;

weff=inf;
k=pi/weff;
Ky=k;

damping=2e-4;
dH0=0.2*796;

i1=1;

dkx=1e-4;

k_max=0.03;

limitation=0.74;

for kx=0:dkx:k_max
    
    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f3=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h)).*exp(i*ky*d)./(2*pi);
    f4=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h)).*exp(i*ky*d)./(2*pi);

    
    Fkxyy0(i1)=integral(f1,-limitation,limitation);
    Fkxzz0(i1)=integral(f2,-limitation,limitation);
    Fkxyyd(i1)=real(integral(f3,-limitation,limitation));
    Fkxzzd(i1)=real(integral(f4,-limitation,limitation));

    Omiga_yy(i1)=(Wh+Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1)));
    Omiga_zz(i1)=(Wh+Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1)));
    
    W0(i1)=sqrt(Omiga_yy(i1).*Omiga_zz(i1));
    
    W1(i1)=sqrt((Omiga_yy(i1)+Wm.*Fkxyyd(i1)).*(Omiga_zz(i1)+Wm.*Fkxzzd(i1)));
    W2(i1)=sqrt((Omiga_yy(i1)-Wm.*Fkxyyd(i1)).*(Omiga_zz(i1)-Wm.*Fkxzzd(i1)));
    
   
    dwm_dwh(i1)=(Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+2.*Wh)./(2.*W0(i1));
    
    i1=i1+1;
end


k1=0:dkx:k_max;

ff0=W0'./(2*pi);
ff1=W1'./(2*pi);
ff2=W2'./(2*pi);

figure
plot(k1,ff0,'b',k1,ff1,'r',k1,ff2,'k')
%axis([-0.030 0.03 0 20])




