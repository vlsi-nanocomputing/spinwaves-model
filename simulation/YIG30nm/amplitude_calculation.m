clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;                                    %thinckness (nm)
w=30;                                    %width(nm)
SW_frequency=2.29; %  [GHz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dkx=1e-4;
kmax=0.025;
kmin=0.001;


limitation=10;
i1=1;
for kx=dkx:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limitation,limitation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limitation,limitation);                                                                                                       %integral

    
    wm0(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));

    
    % nonlinear shift coefficient Tkx calculation
    Akx(i1) = Wh+Wm/2.*(2*(Le.^2).*(kx.^2+Ky.^2)+Fkxyy0(i1)+Fkxzz0(i1));  % [rad/ns]
    Bkx(i1) = Wm/2.*(Fkxyy0(i1)-Fkxzz0(i1));  % [rad/ns]
    

    i1=i1+1;
end
    

k1 = dkx:dkx:kmax;
ff0 = wm0/(2*pi);

%%% amplitude ak   
ukx = sqrt((Akx+wm0)./(2*wm0));
vkx = -sign(Bkx).*sqrt((Akx-wm0)./(2*wm0));

i1=1;
N=size(dkx:dkx:kmax);
N=N(2);
for i1=1:1:N  % for each kx 
    akx = linspace(0,0.8);
    Mz = 2400;
    Mz_Ms_0 = Mz/Ms;  % Mz/Ms
    Mz_Ms = akx.*sqrt(2-abs(akx).^2).*(ukx(i1)-vkx(i1));
    akx_vec(i1) = interp1(Mz_Ms,akx,Mz_Ms_0); % akx(kx) amplitude as a function of kx
end


interp1(k1,akx_vec,interp1(ff0,k1,SW_frequency)) % 2.39 GHz

figure
plot(k1,akx_vec)
ylabel('amplitude(kx)')

% akx=0.0138 is the SW amplitude excited by 4mT. This is the input of DC2
% when at the DC1 input we have A='1' and B='1'. In this case, at the DC2 input
% there is a SW of 2*a_SW/sqrt(2) (constructive interference), where a_SW is the
% amplitude of logic '1'.
%        0.0138 = 2*a_SW/sqrt(2)
%          a_SW = 0.00976
