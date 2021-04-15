clear all
close all
clc

x=linspace(0,100,1000);
Lc=12; % um
x_freepath = 88; %um

P1 = cos(x*pi/(2*Lc)).^2;
P2 = sin(x*pi/(2*Lc)).^2;

P1_real = cos(x*pi/(2*Lc)).^2.*exp(-2*x/x_freepath);
P2_real = sin(x*pi/(2*Lc)).^2.*exp(-2*x/x_freepath);



figure
[ax1]=plotyy(x,P1,x,P2);
xlabel('Distance x  [um]','FontSize',20)
ylabel(ax1(1),'Normalized power in the first waveguide','FontSize',15)
ylabel(ax1(2),'Normalized power in the second waveguide','FontSize',15)
title('Power transfer in dipolarly coupled waveguides','FontSize',15)
grid on


figure
hold on
plot(x,P1_real)
plot(x,P2_real)
hold off
grid on
xlabel('Distance x  [um]','FontSize',20)
legend('Normalized power in the first waveguide','Normalized power in the second waveguide')
title('Power transfer in dipolarly coupled waveguides','FontSize',15)
