function [sum_bit,carry_bit] = HA_ver1(in_A,in_B)
% Half-adder


[DC1_out,DC1_out_I] = DC1_ver1(in_A,in_B);
[sum_bit,carry_bit] = DC2_ver1(DC1_out);

sum_bit = amplifier_ver1(sum_bit,4); 

end

