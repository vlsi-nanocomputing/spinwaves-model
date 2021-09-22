function [bin_number] = dec_to_bin(dec_number, Nbit)
% the function convert a decimal number to unsigned representation
if dec_number >= 2^Nbit
    disp('dec_to_bin error: the number of bits is not enough to represent the required number')
else
    n = 1;
    for j = (Nbit-1):-1:0
        if dec_number >= 2^j
            bin_number(n) = 1;
            dec_number = dec_number - 2^j;
        else
            bin_number(n) = 0;
        end
        n = n+1;
    end
end


end

