function [] = signal_plot(in_signal)

t = linspace(0,0.00001,1000);
y = in_signal(1) * sin(2*pi*in_signal(2)*t*10^9 + in_signal(3));

plot(t,y)
end

