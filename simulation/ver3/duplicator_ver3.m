function [out_signal1,out_signal2] = duplicator_ver3(in_signal)
% The function receives an input signal, and it gives two same signals at
% the output.(The 3 signals are same)

cd common
SW_parameters % script
cd ..

in_signal2 = phase_shifter_ver3(DAC_ver3(0),pi/2);
gain1 = 1/0.4924;
gain2 = 1/(1-0.4924);

[out_signal1, out_signal2] = DC1_ver3(in_signal, in_signal2);

out_signal1 = amplifier_ver3(out_signal1,gain1);
out_signal2 = amplifier_ver3(out_signal2,gain2);
out_signal2 = phase_shifter_ver3(out_signal2,pi/2);

end

