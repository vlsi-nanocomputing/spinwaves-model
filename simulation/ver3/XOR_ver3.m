function [XOR_out] = XOR_ver3(in_A,in_B)


in_B = phase_shifter_ver3(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1_ver3(in_A,in_B);
[out_S,out_C] = DC2_ver3(DC1_out);


XOR_out = regenerator_S_ver3(out_S);

end

