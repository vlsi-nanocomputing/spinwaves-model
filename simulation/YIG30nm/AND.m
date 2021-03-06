function [AND_out] = AND(in_A,in_B)
in_B = phase_shifter(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1(in_A,in_B);
[out_S, AND_out] = DC2(DC1_out);

% % 
AND_out = regenerator_C(AND_out);

end
