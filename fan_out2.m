function [out1,out2] = fan_out2(in_signal)

[out1,out2] = DC1(in_signal,0);

out1 = amplifier(out1);
out2 = amplifier(out2);

end

