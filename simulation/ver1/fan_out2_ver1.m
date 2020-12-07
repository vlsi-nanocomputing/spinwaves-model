function [out1,out2] = fan_out2_ver1(in_signal)

[out1,out2] = DC1_ver1(in_signal,0);

out1 = amplifier_ver1(out1);
out2 = amplifier_ver1(out2);

end

