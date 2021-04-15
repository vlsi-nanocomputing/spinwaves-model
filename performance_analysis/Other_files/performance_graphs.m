clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%% 28nm no load %%%%%%%%%%%%%%%%%%%%%%%%%
% Area_CMOS = [5.22, 10.44, 20.88, 41.78, 83.56];  % [um^2]
% tpd_CMOS = [0.58, 0.66, 0.81, 1.12, 1.73];  % [ns]
% Power_CMOS = [2.11e-4, 4.53e-4, 9.49e-4, 1.957e-3, 3.97e-3];  % [mW]
% Energy_CMOS = Power_CMOS .* tpd_CMOS * 1000;  % [fJ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% 15nm no load %%%%%%%%%%%%%%%%%%%%%%%%%
frequency = 1/10; % GHz,  T_clk = 10ns
Area_CMOS_15nm = [2.36, 4.72, 9.44, 18.87, 37.75];  % [um^2]
tpd_CMOS_15nm = [0.0223, 0.0401, 0.0757, 0.1469, 0.2893];  % [ns]
static_power_CMOS_15nm = [58.3, 116.5, 233.2, 466.2, 932.4];  % [nW]
dynamic_power_CMOS_15nm = [67, 148.2, 315.6, 657.6, 1340.5];  % [nW]
Energy_CMOS_15nm = (static_power_CMOS_15nm .* tpd_CMOS_15nm + dynamic_power_CMOS_15nm/frequency)*1e-3; % [fJ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






addpath('..')
addpath('../YIG')
addpath('../Circuits')
addpath('../YIG/YIG_100nm')
model = 2;
N = 1;
for i=1:5
    N = N*2;
    RCA_nbit
    Area_100nm(i) = area_RCA_nbit;
    Energy_100nm(i) = energy_RCA_nbit;
    tpd_100nm(i) = tpd_RCA_nbit;
end

%%%%%%%%%%%%%%%%%%% YIG 30nm %%%%%%%%%%%%%%%%%%%
addpath('../YIG/YIG_30nm')
model = 3;
N = 1;
for i=1:5
    N = N*2;
    RCA_nbit
    Area_30nm(i) = area_RCA_nbit;
    Energy_30nm(i) = energy_RCA_nbit;
    tpd_30nm(i) = tpd_RCA_nbit;
end







N = [2 4 8 16 32];
figure
hold on
plot(N,Area_100nm,'LineWidth',2)
plot(N,Area_30nm,'LineWidth',2)
plot(N,Area_CMOS_15nm,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Area [um^2]','FontSize',20)
legend('YIG (100nm)','YIG (30nm)', 'CMOS (15nm)')

figure
hold on
plot(N,Energy_100nm*10^(-3),'LineWidth',2)
plot(N,Energy_30nm*10^(-3),'LineWidth',2)
plot(N,Energy_CMOS_15nm,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Energy [fJ/operation]','FontSize',20)
legend('YIG (100nm)','YIG (30nm)', 'CMOS (15nm)')

figure
hold on
plot(N,tpd_100nm,'LineWidth',2)
plot(N,tpd_30nm,'LineWidth',2)
plot(N,tpd_CMOS_15nm,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Propagation delay [ns]','FontSize',20)
legend('YIG (100nm)','YIG (30nm)', 'CMOS (15nm)')

ax = figures(3).CurrentAxes;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;