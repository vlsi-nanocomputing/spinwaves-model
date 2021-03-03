function [XOR_out] = XOR(in_A,in_B)


in_B = phase_shifter(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1(in_A,in_B);
[out_S,out_C] = DC2(DC1_out);


XOR_out = regenerator_S(out_S);

end

