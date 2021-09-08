function [] = signal_plotting(in_signal,model,varargin)

N_signal = size(in_signal);
N_signal = N_signal(1);
t_max = max(in_signal(:,4)) + 10/in_signal(1,2); % ns

SW_parameters

figure
hold on
for i=1:N_signal

    a0 = in_signal(i,1);
    f = in_signal(i,2);         % GHz
    phase = in_signal(i,3);     % rad
    tpd = in_signal(i,4);       % ns
    
  
    t1 = linspace(0,tpd,2);
    y1 = 0*t1;
    t2 = [tpd:0.005:t_max];
    y2 = a0*sin(2*pi*f*(t2-tpd)+phase);

    t = [t1 t2];
    y = [y1 y2]-(i-1)*2.2*SW_amplitude;  % offset to distinguish the different signals

    plot(t,y)
end
hold off
xlabel('time [ns]','FontSize',20)
grid on
axis([0, t_max, -2.2*SW_amplitude*N_signal+1.1*SW_amplitude, 1.1*SW_amplitude])
if nargin > 2  % optional legend
    legend(varargin{1:end})
end

end