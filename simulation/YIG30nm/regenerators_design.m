clear all
close all
clc

cd common
SW_parameters
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 10;                                    %thinckness (nm)
w = 30;                                    %width(nm)
gap = 10;
d = w+gap;
B = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % output S
% gain_in = 8;
% a10 = 0.042919800099790 * sqrt(gain_in) ;
% a01 = 0.045715737312973 * sqrt(gain_in) ;
% a11 = 0.013812976999434 * sqrt(gain_in) ;

% % output S: after 900nm
% gain_interm = 2;
% a10 = 0.084674257671247 * sqrt(gain_interm) ;
% a01 = 0.096033838031988 * sqrt(gain_interm) ;
% a11 = 0.000040400167428 * sqrt(gain_interm) ;

% output C
a10 = 0.015826760179143;
a01 = 0.013545867914880;
a11 = 0.092326830728056;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
delta_ph10=0;
delta_ph01=0;
delta_ph11=0;

DC_design = [h, w, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation2, DC_design);
cd ..

ff1 = wm1/(2*pi);
ff2 = wm2/(2*pi);

Lw_max = 1450;
Lw=1:1:Lw_max;
dl = 1;
for i1=1:dl:Lw_max
    a10 = a10*exp(-dl/x_freepath);
    a01 = a01*exp(-dl/x_freepath);
    a11 = a11*exp(-dl/x_freepath);

    ff1_s10 = ff1+Tkx.*abs(a10).^2;
    ff2_s10 = ff2+Tkx.*abs(a10).^2;
    
    ff1_s01 = ff1+Tkx.*abs(a01).^2;
    ff2_s01 = ff2+Tkx.*abs(a01).^2;
    
    ff1_s11 = ff1+Tkx.*abs(a11).^2;
    ff2_s11 = ff2+Tkx.*abs(a11).^2;
    
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


out10 = a10*sqrt((pow_par10(end)));
out01 = a01*sqrt((pow_par01(end)));
out11 = a11*sqrt((pow_par11(end)));
% normalization([out10,out01,out11])

% x = linspace(0,1,1000);
% y = x*(a10/a01)^2;
% figure
% hold on
% plot(x,y)
% plot(1-x,1-y)
% hold off
% xlabel('01','FontSize',20)
% ylabel('10','FontSize',20)
% legend('waveguide1 output','waveguide2 output')

figure
hold on
plot(Lw,pow_par10,'LineWidth',1.5)
plot(Lw,pow_par01,'LineWidth',1.5)
plot(Lw,pow_par11,'LineWidth',1.5)
hold off
xlabel('L_w  [nm]','FontSize',20)
legend('10','01','11')
% legend('10','01')
title('Normalized output power','FontSize',15)
