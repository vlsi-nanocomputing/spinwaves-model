function [out,out_I] = DC1_ver1(in_A,in_B)
% Directional Coupler 1

% in_A and in_B are the input spin waves of the circuit
% DC1 works in a linear regime and acts as a symmetric power splitter for
% each of the two inputs. out=(in_A+in_B)/2

out = (in_A+in_B)/2;

% the out_I is an additional output, that can be used to create a second "out"
out_I = (in_A+in_B)/2;

end

