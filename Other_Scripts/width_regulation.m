clear
close all
clc
%%%% input:     setting section
%%%% output:    width of waveguide, group velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setting section %%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;                % thinckness (nm)
w=100;                % initial width(nm)
delta_k = -7*9.9331e-05;  % (rad/nm), variation with respect to the k(w): positive or negative
SW_frequency = 2.282; % GHz
width_min = 80;      % nm
width_max = 110;      % nm
% the program searches the width value in the range from width_min to width_max
resolution=1;        % resolution to find new width [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%w=25;  % k=9.69, vg=156

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% physical parameters %%%%%%%%%%%%%%%%%%%%
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


dkx=1e-5;
kmax=0.02;
limiation=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initial k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i1=1;
for kx=dkx:dkx:kmax                       

    f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
    f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

    Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
    Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral

    wm(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));

    i1=i1+1;
end
k1=dkx:dkx:kmax;
ff=wm./(2*pi);
init_k = interp1(ff,k1,SW_frequency);
target_k = init_k + delta_k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


init_w = w;
err_new=1000;
vg_new=0;
flag = 0;
for w = width_min : resolution : width_max
    i1=1;
    for kx=dkx:dkx:kmax                       

        f1=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*ky.^2./(kx.^2+ky.^2).*(1-(1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);
        f2=@(ky)(abs(2*sqrt(2./(1+sinc(k*w./pi))).*(ky.*cos(k.*w./2).*sin(ky.*w./2)-k.*(cos(ky.*w./2).*sin(k.*w./2)))./(ky.^2-k.^2))).^2./w.*((1-exp(-sqrt(kx.^2+ky.^2).*h))./(sqrt(kx.^2+ky.^2).*h))./(2*pi);

        Fkxyy0(i1)=integral(f1,-limiation,limiation);                                                                                                       %integral
        Fkxzz0(i1)=integral(f2,-limiation,limiation);                                                                                                       %integral


        wm(i1)=sqrt((Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))).*(Wh+Wm*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))));
        %dwm_dwh(i1)=(Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxyy0(i1))+Wm.*(Le.^2.*(kx.^2+Ky.^2)+Fkxzz0(i1))+2.*Wh)./(2.*wm(i1));

        i1=i1+1;

    end


     k1=dkx:dkx:kmax;
     ff=wm./(2*pi);

     dw=diff(wm);
     dk=diff(k1);
     N_length=length(dk);

     velocity=abs(dw./dk/1000);
     %lifetime=1./((damping+r*dH0*1e-9./(2*wm)).*wm.*dwm_dwh)';
     %distance=abs(lifetime(1:N_length).*velocity');

%%%%%%%%%%%     
     err_old = err_new;
     err_new = abs(interp1(ff,k1,SW_frequency) - target_k);
     vg=vg_new;
     vg_new = interp1(k1(1:N_length),velocity,interp1(ff,k1,SW_frequency));
     if err_new>err_old
        w=w-resolution;
        fprintf('Starting from the waveguide with w = %d nm, to obtain a k variation delta_k = %d rad/nm you can use a waveguide with w = %d nm. \n',init_w,delta_k,w)
        fprintf('Using this new width you can obtain the following parameters: wavenumber k = %d rad/nm, and group velocity vg = %d km/s',err_old + target_k, vg)
        flag = 1;
        break
     end
end

if flag == 0
    fprintf('Starting from the waveguide with w = %d nm, it is not possible to find a width to have a delta_k = %d rad/nm in the range [%dnm, %dnm] \n', init_w,delta_k,width_min,width_max)
end

