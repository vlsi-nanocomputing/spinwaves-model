function [sum_bit,carry_bit] = HA(in_A,in_B)
% Half-adder


[DC1_out,DC1_out_I] = DC1(in_A,in_B);
[sum_bit,carry_bit] = DC2(DC1_out);

sum_bit = amplifier(sum_bit); 

end

