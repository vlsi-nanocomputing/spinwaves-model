function [out_S,out_C] = HA_ver3(in_A,in_B)

in_B = phase_shifter_ver3(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1_ver3(in_A,in_B);
[out_S,out_C] = DC2_ver3(DC1_out);

% % 
out_S = regenerator_S_ver3(out_S);
out_C = regenerator_C_ver3(out_C);

end

