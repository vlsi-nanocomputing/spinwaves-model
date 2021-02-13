function [AND_out] = AND_ver3(in_A,in_B)
in_B = phase_shifter_ver3(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1_ver3(in_A,in_B);
[out_S, AND_out] = DC2_ver3(DC1_out);

% % 
AND_out = regenerator_C_ver3(AND_out);

end
