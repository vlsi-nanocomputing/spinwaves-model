function [X] = CLA_4bit_ver2(A,B,carry)


X = zeros(5,3);
p_i = zeros(4,3);
gi = zeros(4,3);
ci = zeros(4,3);
ai_xor_bi = zeros(4,3);



for i=1:4
	[gi(i,:), p_i(i,:), ai_xor_bi(i,:)] = CLA_block_A1_ver2(A(i,:),B(i,:));
end
[g0_1, g0_2] = duplicator_ver2(gi(4,:));
[p0_1, p0_2] = duplicator_ver2(p_i(4,:));

% block B1 for a0a1-b0b1
[G01, P01] = CLA_block_B1_ver2(gi(3,:), p_i(3,:), g0_1, p0_1);

% final block B2 
[ci(4,:), ci(2,:)] = CLA_block_B2_ver2(carry, G01, P01);

% B2 of the second level
[ci(4,:), ci(3,:)] = CLA_block_B2_ver2(ci(4,:), g0_2, p0_2);
[ci(2,:), ci(1,:)] = CLA_block_B2_ver2(ci(2,:), gi(2,:), p_i(2,:));


[ci(1,:), c3] = duplicator_ver2(ci(1,:));
% A2 blocks to generate the sum bits
for i=1:4
	X(i+1,:) = CLA_block_A2_ver2(ai_xor_bi(i,:), ci(i,:));
end

X(1,:) = OR_ver2( AND_ver2(c3,p_i(1,:)), gi(1,:) );


end

