% Generic simulator
clear all 

%%%%%%%%%%%%%%%%%%%%%%%% simulation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_bin = [0 1;1 1]; % binary representation of the input A(from MSB to LSB)
B_bin = [1 0;0 1]; % binary representation of the input B(from MSB to LSB)
C_bin = [1;0]; % binary representation of the input C(from MSB to LSB)
D_bin = [1 1;0 1]; % binary representation of the input D(from MSB to LSB)

Nbit = 2;    % it is used by RCA and CSA (parallelism). For the Carry-Skip Adder, the Nbit must be a multiple of 4, which is a constraint of the CSA model
opt_parameters = {'out_signal_plot'}; % optional parameters, it can be empty
titleFontSize = 25;   % title FontSize of the plots
axisFontSize = 13;    % axes FontSize of the plots
labelFontSize = 20;   % labels FontSize of the plots
legendFontSize = 24;  % legend FontSize of the plots
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

%%%%%%%%%%%%%%%%%%%%% Choosing of the simulation circuit %%%%%%%%%%%%%%%%
circuit=0;
while circuit ~= [1,2,3,4,5,6,7,8,9,10]
    fprintf('\nChoose one simulation circuit from the following list:');
    fprintf('\n  1) AND(A,B)');
    fprintf('\n  2) AND4(A,B,C,D)');
    fprintf('\n  3) %d-bit Carry_skip_adder(A,B,C=Carry_in)',Nbit);
    fprintf('\n  4) FA(A,B,C)');
    fprintf('\n  5) HA(A,B)');
    fprintf('\n  6) Mux2to1(A,B,C=sel)');
    fprintf('\n  7) NOT(A)');
    fprintf('\n  8) OR(A,B)');
    fprintf('\n  9) %d-bit RCA(A,B,C=Carry_in)',Nbit);
    circuit = input('\n  10) XOR(A,B) \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_A = size(A_bin);
N_simulation = size_A(1);
for ii = 1:N_simulation
    switch circuit
        case 1
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            AND_out = AND(A, B, model, opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 2
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            C = DAC(C_bin(ii,:),model);
            D = DAC(D_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('", C = "')
            fprintf('%d',C_bin(ii,:))
            fprintf('", D = "')
            fprintf('%d',D_bin(ii,:))
            fprintf('"\n')
            AND4_out = AND4(A, B, C, D, model, opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 3
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            C = DAC(C_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('", C = "')
            fprintf('%d',C_bin(ii,:))
            fprintf('"\n')
            CSA_out = carry_skip_adder(A,B,C,Nbit,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 4
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            C = DAC(C_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('", C = "')
            fprintf('%d',C_bin(ii,:))
            fprintf('"\n')
            [FA_S,FA_C] = FA(A,B,C,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 5
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            [HA_S,HA_C] = HA(A, B, model, opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 6
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            C = DAC(C_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('", C = "')
            fprintf('%d',C_bin(ii,:))
            fprintf('"\n')
            MUX2TO1_out = mux2to1(A,B,C,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 7
            A = DAC(A_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('"\n')
            NOT_out = NOT(A,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 8
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            OR_out = OR(A,B,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 9
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            C = DAC(C_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('", C = "')
            fprintf('%d',C_bin(ii,:))
            fprintf('"\n')
            RCA_out = RCA_Nbit(A,B,C,Nbit,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
        case 10
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            XOR_out = XOR(A,B,model,opt_parameters{:})
            title('Simulation ' + string(ii), 'FontSize', titleFontSize)
            ax = gca;
            ax.Legend.FontSize = legendFontSize;
            ax.XAxis.FontSize = axisFontSize;
            ax.YAxis.FontSize = axisFontSize;
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
    end
end



end % if model~=4