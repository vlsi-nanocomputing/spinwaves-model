% simulation
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%% digital input setting %%%%%%%%%%%%%%%%%%%%%%
A = 1;
B = 0;

C = 1;


% OR_ver1(1,1) % fast verification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% digital to analog conversion %%%%%%%%%%%%%%%%%%%%%%%
err = 0;

switch A
    case 0
        in_A(1) = 0;
    case 1
        in_A(1) = 0.153;
    otherwise
        display("Error: check your digital inputs")
        err = 1;
end
in_A(2) = 2.282;
in_A(3) = 0;

switch B
    case 0
        in_B(1) = 0;
    case 1
        in_B(1) = 0.153;
    otherwise
        display("Error: check your digital inputs")
        err = 1;
end
in_B(2) = 2.282;
in_B(3) = 0;

switch C
    case 0
        in_C(1) = 0;
    case 1
        in_C(1) = 0.153;
    otherwise
        display("Error: check your digital inputs")
        err = 1;
end
in_C(2) = 2.282;
in_C(3) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation circuit %%%%%%%%%%%%%%%%%%%%%%%%%
% % if err == 0
% %    OR_out = OR_ver2(in_A,in_B);
% %    OR_out = OR_ver2(OR_out,in_B);
% % end
% % 
% % OR_out(1)^2/0.153^2

% % NOT_ver2(in_A)

% % mux2to1_ver2(in_A,in_B,in_C)

XOR_ver2(in_A,in_B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


