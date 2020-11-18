
clear all
clc

A1=[0.153, 2.282, 0];
B1=[0, 2.282, 0];
[C1,D1] = HA_without_amp_ver2(A1,B1);

A2=[0, 2.282, 0];
B2=[0.153, 2.282, 0];
[C2,D2] = HA_without_amp_ver2(A2,B2);


A3=[0.153, 2.282, 0];
B3=[0.153, 2.282, 0];
[C3,D3] = HA_without_amp_ver2(A3,B3);

i1=1;
for gain = 1.8:0.01:2.2
    err_sum1(i1) = C1(1).^2*gain/0.153.^2;
    err_sum2(i1) = C2(1).^2*gain/0.153.^2;
    err_sum3(i1) = C3(1).^2*gain/0.153.^2;
    i1=i1+1;
end


gain = [1.8:0.01:2.2];
points = size(gain);
points = points(2);
hold on
plot(gain,err_sum1,'LineWidth',1.5);
plot(gain,err_sum2,'LineWidth',1.5);
plot(gain,err_sum3,'LineWidth',1.5);
plot(gain,ones(1,points),'LineWidth',1.5)
hold off
axis([1.8 2.2 0 1.1])
xlabel('gain','FontSize',20)
legend('err\_sum1','err\_sum2','err\_sum3','100%')










% 
% i1=1;
% for gain = 1.8:0.01:2.2
%     err_sum1(i1) = C1(1).^2*gain-0.153.^2;
%     err_sum2(i1) = C2(1).^2*gain-0.153.^2;
%     err_sum3(i1) = C3(1).^2*gain-0;
%     i1=i1+1;
% end
% 
% gain = [1.8:0.01:2.2];
% points = size(gain);
% points = points(2);
% hold on
% plot(gain,err_sum1,'LineWidth',1.5);
% plot(gain,err_sum2,'LineWidth',1.5);
% plot(gain,err_sum3,'LineWidth',1.5);
% plot(gain,0*ones(1,points),'LineWidth',1.5)
% hold off
% axis([1.8 2.2 -0.01 0.01])
% xlabel('gain','FontSize',20)
% legend('err\_sum1','err\_sum2','err\_sum3','SW amplitude','0')