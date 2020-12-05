function [] = signal_plot(in_signal)

N_signal = size(in_signal);
N_signal = N_signal(1);
t_max = max(in_signal(:,4)) + 3/in_signal(1,2); % ns


hold on
for i=1:N_signal

    a0 = in_signal(i,1);
    f = in_signal(i,2);         % GHz
    phase = in_signal(i,3);     % rad
    tpd = in_signal(i,4);       % ns
    
  
    t1 = linspace(0,tpd,2);
    y1 = 0*t1;
    t2 = linspace(tpd,t_max,1000);
    y2 = a0*sin(2*pi*f*(t2-tpd)+phase);

    t = [t1 t2];
    y = [y1 y2];

    plot(t,y)
    
    
end
hold off
xlabel('time [ns]','FontSize',20)


end