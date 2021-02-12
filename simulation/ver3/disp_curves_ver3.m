clear all
close all
clc



%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;         % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
DC2_akx=0.0138/sqrt(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;

figure
hold on
for i1=0:1:1

    limitation=0.74-i1*0.4;

    DC2_design = [h, w, d, B];
    cd common
    SW_parameters
    [wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
    cd ..

    k1=dkx:dkx:kmax;
    DC2_ff1=wm1./(2*pi);
    DC2_ff2=wm2./(2*pi);
    
    plot(k1,DC2_ff1,'LineWidth',2)
    plot(k1,DC2_ff2,'LineWidth',2)
end
p=size(k1);
p=p(2);
plot(k1,2.282*ones(1,p),'LineWidth',2)
hold off
legend('f1\_0.74','f2\_0.74','f1\_10','f2\_10','2.282 GHz')



