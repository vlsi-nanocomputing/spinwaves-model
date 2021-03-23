% Generic simulator
clear all 

%%%%%%%%%%%%%%%%%%%%%%%% simulation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_bin = [1]; % binary representation of the input A(from MSB to LSB)
B_bin = [1]; % binary representation of the input B(from MSB to LSB)
C_bin = [0]; % binary representation of the input C(from MSB to LSB)
D_bin = [1]; % binary representation of the input D(from MSB to LSB)

Nbit = 4;    % it is used by RCA and CSA (parallelism). For the CSA, the Nbit must be a multiple of 4, that is a constraint of the CSA model
opt_parameters = {'regS_dispersion_curves','DC2_Lc_avg','DC1_dispersion_curves','regC_Lc_avg'}; % optional parameters, it can be empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Choosing of model category %%%%%%%%%%%%%%%%%%%
model=0;
while model~=1 && model~=2 && model~=3 && model~=4 
model = input('Choose one model category from the following list:\n 1) YIG (100 nm) Behavioral Model \n 2) YIG (100 nm) Physical Model \n 3) YIG (30 nm) Physical Model \n 4) QUIT \n');
end
if model~=4

switch model 
    case 1
        model_path = 'Building_blocks/YIG100nm_Behavioral_model';
    case 2
        model_path = 'Building_blocks/YIG100nm_Physical_model';
    case 3
        model_path = 'Building_blocks/YIG30nm_Physical_model';
end
addpath(model_path)
addpath('Building_blocks/Common')
addpath('Circuits')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% digital to analog conversion %%%%%%%%%%%%%%%%%%%%
A = DAC(A_bin,model);
B = DAC(B_bin,model);
C = DAC(C_bin,model);
D = DAC(D_bin,model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Choosing of the simulation circuit %%%%%%%%%%%%%%%%
circuit=0;
while circuit ~= [1,2,3,4,5,6,7,8,9,10]
    fprintf('\nChoose one simulation circuit from the following list:');
    fprintf('\n  1) AND(A,B)');
    fprintf('\n  2) AND4(A,B,C,D)');
    fprintf('\n  3) %d-bit Carry_skip_adder(A,B,Carry_in)',Nbit);
    fprintf('\n  4) FA(A,B,C)');
    fprintf('\n  5) HA(A,B)');
    fprintf('\n  6) Mux2to1(A,B,sel)');
    fprintf('\n  7) NOT(A)');
    fprintf('\n  8) OR(A,B)');
    fprintf('\n  9) %d-bit RCA(A,B,Carry_in)',Nbit);
    circuit = input('\n  10) XOR(A,B) \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch circuit
    case 1
        AND_out = AND(A, B, model, opt_parameters{:})
    case 2
        AND4_out = AND4(A, B, C, D, model, opt_parameters{:})
    case 3
        CSA_out = carry_skip_adder(A,B,C,Nbit,model,opt_parameters{:})
    case 4
        [FA_S,FA_C] = FA(A,B,C,model,opt_parameters{:})
    case 5
        [HA_S,HA_C] = HA(A, B, model, opt_parameters{:})
    case 6
        MUX2TO1_out = mux2to1(A,B,C,model,opt_parameters{:})
    case 7
        NOT_out = NOT(A,model,opt_parameters{:})
    case 8
        OR_out = OR(A,B,model,opt_parameters{:})
    case 9
        RCA_out = RCA_Nbit(A,B,C,Nbit,model,opt_parameters{:})
    case 10
        XOR_out = XOR(A,B,model,opt_parameters{:})
end

end % if model~=4