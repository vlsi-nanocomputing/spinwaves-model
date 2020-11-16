function [out_S,out_C] = HA_ver2(in_A,in_B)


[DC1_out,DC1_out_I] = DC1_ver2(in_A,in_B);
[out_S,out_C] = DC2_ver2(DC1_out);

out_S(1) = out_S(1); 
out_C(1) = out_C(1)/sqrt(2); % damping  


end

