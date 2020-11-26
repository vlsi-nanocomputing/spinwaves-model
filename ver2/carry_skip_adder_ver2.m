function [X] = carry_skip_adder_ver2(A,B,carry,Nbit)
% N-bit carry skip adder
% Please use a N as a multiple of 4
% NB: the function describes the behavioral simulation of the device 
% (carry skip adder), the following structure is different with respect to
% the VHDL structure.

N_block = Nbit/4;

%%% We need to duplicate A and B, where the first ones are sent to the
% FAs, and second ones are used to calculate the propagates.
A1 = zeros(Nbit,3);
A2 = zeros(Nbit,3);
B1 = zeros(Nbit,3);
B2 = zeros(Nbit,3);
for i=1:Nbit
    [A1(i,:),A2(i,:)] = duplicator_ver2(A(i,:));
    [B1(i,:),B2(i,:)] = duplicator_ver2(B(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% propagate calculation %%%%%%%%%%%%%%%%%%%%%%%%%
p = zeros(Nbit,3); % propagates of every (ai,bi)
for i=1:Nbit
    p(i,:) = XOR_ver2( A1(i,:), B1(i,:) );
end
P_block = zeros(N_block,3); % propagate of every block
for i=1:N_block
    P_block(i,:) = AND4_ver2( p((i-1)*4+1,:), p((i-1)*4+2,:), p((i-1)*4+3,:), p((i-1)*4+4,:) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% 4-bit RCA of every block %%%%%%%%%%%%%%%%%%%%%
if rem(Nbit,4) == 0
    for i=N_block:-1:1 % for every block
        RCA_out = RCA_Nbit_ver2( A2((i-1)*4+1:(i-1)*4+4, :), B2((i-1)*4+1:(i-1)*4+4, :), carry, 4 );
        X( (i-1)*4+2:(i-1)*4+5, :) = RCA_out(2:end, :);
        carry = mux2to1_ver2(RCA_out(1,:), carry, P_block(i,:));
    end
    X(1,:) = carry;
else
    display('please check your Nbit')
end

