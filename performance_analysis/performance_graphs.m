clear all
close all
clc

Area_CMOS = [5.22, 10.44, 20.88, 41.78, 83.56];  % [um^2]
tpd_CMOS = [0.58, 0.66, 0.81, 1.12, 1.73];  % [ns]
Power_CMOS = [2.11e-4, 4.53e-4, 9.49e-4, 1.957e-3, 3.97e-3];  % [mW]
Energy_CMOS = Power_CMOS .* tpd_CMOS * 1000;  % [fJ]

Area = zeros(1,5);
Energy = zeros(1,5);
tpd = zeros(1,5);
N = 1;
for i=1:5
    N = N*2;
    RCA_nbit
    Area(i) = A_RCA_nbit;
    Energy(i) = E_RCA_nbit;
    tpd(i) = tpd_RCA_nbit;
end
N = [2 4 8 16 32];

figure
hold on
plot(N,Area,'LineWidth',2)
plot(N,Area_CMOS,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Area [um^2]','FontSize',20)
legend('YIG (100nm)', 'CMOS (28nm)')

figure
hold on
plot(N,Energy*10^(-3),'LineWidth',2)
plot(N,Energy_CMOS,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Energy [fJ/operation]','FontSize',20)
legend('YIG (100nm)', 'CMOS (28nm)')

figure
hold on
plot(N,tpd,'LineWidth',2)
plot(N,tpd_CMOS,'LineWidth',2)
hold off
xlabel('N-bit','FontSize',20)
ylabel('Propagation delay [ns]','FontSize',20)
legend('YIG (100nm)', 'CMOS (28nm)')


