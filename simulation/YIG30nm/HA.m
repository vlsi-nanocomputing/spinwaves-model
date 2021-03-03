function [out_S,out_C] = HA(in_A,in_B)

in_B = phase_shifter(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1(in_A,in_B);
[out_S,out_C] = DC2(DC1_out);

% % 
out_S = regenerator_S(out_S);
out_C = regenerator_C(out_C);

end
