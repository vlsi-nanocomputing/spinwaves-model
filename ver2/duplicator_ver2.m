function [out_signal1,out_signal2] = duplicator_ver2(in_signal)
% The function receives an input signal, and it gives two same signals at
% the output.(The 3 signals are same)

in_signal2 = [0, 2.282, pi/2];
gain1 = 1/0.500249501323049;
gain2 = 1/(1-0.500249501323049);

[out_signal1, out_signal2] = DC1_ver2(in_signal, in_signal2);

out_signal1 = amplifier_ver2(out_signal1,gain1);
out_signal2 = amplifier_ver2(out_signal2,gain2);
out_signal2 = phase_shifter_ver2(out_signal2,pi/2);

end

