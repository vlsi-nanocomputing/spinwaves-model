function [out_signal1,out_signal2] = duplicator(in_signal)
% The function receives an input signal, and it gives two same signals at
% the output.(The 3 signals are same)

cd common
SW_parameters % script
cd ..

in_signal2 = phase_shifter(DAC(0),pi/2);
gain1 = 2.25;
gain2 = 2.07;

[out_signal1, out_signal2] = DC1(in_signal, in_signal2);

out_signal1 = amplifier(out_signal1,gain1);
out_signal2 = amplifier(out_signal2,gain2);
out_signal2 = phase_shifter(out_signal2,pi/2);

end

