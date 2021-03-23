function [out_signal1,out_signal2] = duplicator(in_signal,model)
% The function receives an input signal, and it gives two same signals at
% the output.(The 3 signals are same)

SW_parameters % script

in_signal2 = phase_shifter(DAC(0,model),pi/2);


[out_signal1, out_signal2] = DC1(in_signal, in_signal2, model);

out_signal1 = amplifier(out_signal1,duplicator_gain1);
out_signal2 = amplifier(out_signal2,duplicator_gain2);
out_signal2 = phase_shifter(out_signal2,pi/2);

end

