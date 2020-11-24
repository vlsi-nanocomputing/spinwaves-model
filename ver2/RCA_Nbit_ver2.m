function [X] = RCA_Nbit_ver2(A,B,carry,Nbit)
% N-bit ripple carry adder
X = zeros(Nbit+1,3);

for j=0:Nbit-1
    [X(Nbit-j+1,:),carry] = FA_ver2( A(Nbit-j,:), B(Nbit-j,:), carry);
end
X(1,:) = carry; % MSB or final carry_out of the RCA

end