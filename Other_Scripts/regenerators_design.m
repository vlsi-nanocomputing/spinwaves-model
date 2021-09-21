clear all
close all
clc
% This script is used to find the length (Lw) of a signle regenerator DC (S or C).
% At the beginning, you need to set your degrated values. The goal is to
% attenuate (using the nonlinearity of DC) the degrated '0' towards zero,
% but keep the level of the '1'. At the end we can use a small amplifier
% to amplify the '1'.

titleFontSize = 60;   % title FontSize of the plots
axisFontSize = 60;    % axes FontSize of the plots
labelFontSize = 65;   % labels FontSize of the plots
legendFontSize = 65;  % legend FontSize of the plots
line_width = 5;       % LineWidth of the lines 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setting section %%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_in = 1/1.3;%9/4.2; % input amplifier gain
model = 'YIG 100nm'; % =1 for YIG100nm Behavioral model, =2 for YIG100nm Physical model, =3 for YIG30nm Physical model
length_max = 1516; % max length of the DC [nm]
resolution = 1; % discretization resolution [nm]
% amplitude of the output S or C for every input combination (10,01,11)
% A_10 = 6.228145e-02;%0.031993494723368;
% A_01 = 6.045905e-02;%0.032549676365182;
% A_11 = 3.007010e-05;%0.009786218005401;
A_10 = 2.130217e-03;%0.031993494723368;  2.416911e-04;%
A_01 = 2.845638e-04;%0.032549676365182; 8.832202e-04;%
A_11 = 8.057632e-02;%0.009786218005401;   7.779518e-02;%
% A_10 = 7.191374e-02;%0.031993494723368;  2.416911e-04;%
% A_01 = 8.000682e-02;%0.032549676365182; 8.832202e-04;%
% A_11 = 2.197420e-05;%0.009786218005401;   7.779518e-02;%

% input amplifier
A_10 = A_10 * sqrt(gain_in);
A_01 = A_01 * sqrt(gain_in);
A_11 = A_11 * sqrt(gain_in);

% initial parameters of this DC
h=30;                                    %thinckness (nm)
w=100;                                    %width(nm)
gap=10; % the gap between the second coupled waveguides [nm]
B=0; % external field [mT]
% the folder that contains "SW_parameters.m" and "DC_equations.m"
addpath('../simulation/Building_blocks/Common')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%
d=w+gap;
SW_parameters
DC_design = [h, w, d, B];
limitation = model_parameters.limitation2;
[wm1, wm2, Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC_design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% design operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=model_parameters.dkx:model_parameters.dkx:model_parameters.kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);


delta_ph10=0;
delta_ph01=0;
delta_ph11=0;
dl=resolution;
Lw=1:1:length_max;

for i1=1:1:length_max
    A_10 = A_10*exp(-dl/model_parameters.x_freepath);
    A_01 = A_01*exp(-dl/model_parameters.x_freepath);
    A_11 = A_11*exp(-dl/model_parameters.x_freepath);

    ff1_s10 = ff1+Tkx.*abs(A_10).^2;
    ff2_s10 = ff2+Tkx.*abs(A_10).^2;
    
    ff1_s01 = ff1+Tkx.*abs(A_01).^2;
    ff2_s01 = ff2+Tkx.*abs(A_01).^2;
    
    ff1_s11 = ff1+Tkx.*abs(A_11).^2;
    ff2_s11 = ff2+Tkx.*abs(A_11).^2;
    
    ks = interp1(abs(ff1_s10),k1,model_parameters.SW_frequency);
    kas = interp1(abs(ff2_s10),k1,model_parameters.SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph10 = delta_ph10 + delta_k*dl;
    
    ks = interp1(abs(ff1_s01),k1,model_parameters.SW_frequency);
    kas = interp1(abs(ff2_s01),k1,model_parameters.SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph01 = delta_ph01 + delta_k*dl;
    
    ks = interp1(abs(ff1_s11),k1,model_parameters.SW_frequency);
    kas = interp1(abs(ff2_s11),k1,model_parameters.SW_frequency);
    delta_k = abs(ks-kas);
    delta_ph11 = delta_ph11 + delta_k*dl;

    
    pow_par10(i1) = cos(delta_ph10/2)^2;
    pow_par01(i1) = cos(delta_ph01/2)^2;
    pow_par11(i1) = cos(delta_ph11/2)^2;
    
end


out10 = A_10*sqrt((pow_par10(end)));
out01 = A_01*sqrt((pow_par01(end)));
out11 = A_11*sqrt((pow_par11(end)));
fprintf('for the L_DC = length_max, you can obtain the following output signal amplitudes: \n')
fprintf('A=1, B=0: output signal amplitude = %d, the normalized power with respect to the logic ''1'' is %d%% \n',out10,normalization(out10,model_parameters))
fprintf('A=0, B=1: output signal amplitude = %d, the normalized power with respect to the logic ''1'' is %d%% \n',out01,normalization(out01,model_parameters))
fprintf('A=1, B=1: output signal amplitude = %d, the normalized power with respect to the logic ''1'' is %d%% \n',out11,normalization(out11,model_parameters))

vector = abs(A_01*sqrt((pow_par01(:)))-A_10*sqrt((pow_par10(:))));
[m, i] = min(vector)
(0.0779/out10)^2
(0.0779/out01)^2
(0.0779/out11)^2

x = linspace(0,1,1000);
y = x*(A_10/A_01)^2;
figure
plot(x,y)
title('y = x * (A\_10/A\_01)^2')

figure
hold on
plot(Lw,pow_par10,'LineWidth',line_width)
plot(Lw,pow_par01,'LineWidth',line_width)
plot(Lw,pow_par11,'LineWidth',line_width)
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','arial')  % Set it to times
grid on
hold off
xlabel('L_w  [nm]','FontSize',labelFontSize)
lgd = legend('A=1,B=0','A=0,B=1','A=1,B=1');
lgd.FontSize = legendFontSize;
set(gca,'FontSize',axisFontSize)
% clegend('10','01')
%title('Normalized output power (%)','FontSize',titleFontSize)