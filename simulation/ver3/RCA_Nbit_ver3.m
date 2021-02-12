function [X] = RCA_Nbit_ver3(A,B,carry,Nbit)
% N-bit ripple carry adder

cd common
SW_parameters % script
cd ..

X = zeros(Nbit+1,N_inf);

for j=0:Nbit-1
    [X(Nbit-j+1,:),carry] = FA_ver3( A(Nbit-j,:), B(Nbit-j,:), carry);
end
X(1,:) = carry; % MSB or final carry_out of the RCA

end