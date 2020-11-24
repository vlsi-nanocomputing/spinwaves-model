function [X] = RCA_Nbit_ver1(A,B,carry,Nbit)

X = zeros(1,Nbit+1);

for j=0:Nbit-1
    [X(Nbit-j+1),carry] = FA_ver1( A(Nbit-j), B(Nbit-j), carry);
end
X(1) = carry; % MSB or final carry_out of the RCA

end