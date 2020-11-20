function [out_S,out_C] = HA_ver2(in_A,in_B)

in_B = phase_shifter_ver2(in_B,pi/2);

[DC1_out,DC1_out_I] = DC1_ver2(in_A,in_B);
[out_S,out_C] = DC2_ver2(DC1_out);


% out_S = amplifier_ver2(out_S,2.02); 
out_S = regenerator_S_ver2(out_S);
out_C(1) = out_C(1)/sqrt(2); % damping  
out_C = regenerator_C_ver2(out_C);

end

