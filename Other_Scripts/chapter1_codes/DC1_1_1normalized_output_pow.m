clear all
close all
clc

Lw=5;  %um
Lc=linspace(1,4);  %um
norm_p1out=cos(pi*Lw./(2*Lc)).^2;
norm_p2out=sin(pi*Lw./(2*Lc)).^2;


figure
[ax1]=plotyy(Lc,norm_p1out*100,Lc,norm_p2out*100);
xlabel('Coupling length L_{c} [um]','FontSize',15)
ylabel(ax1(1),'P_{1out}/(P_{1out}+P_{2out})','FontSize',15)
ylabel(ax1(2),'P_{2out}/(P_{1out}+P_{2out})','FontSize',15)
title('Normalized output powers (%)','FontSize',20)
grid on



