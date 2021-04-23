% Generic simulator
clear all 

% This interface script allows to simulate all the circuit of "Circuits"
% folder
% To run the script you must set all the parameters of "simulation setting"
% section.


% A_bin, B_bin, C_bin and D_bin are binary inputs of the circuit. You can
% set more than one simulation using multiple rows of the inputs: every row
% is a std_logic_vector(MSB downto LSB) that will be used for a single
% simulation. In that case, the code runs a number of simulations
% equal to the number of the input A (rows).

% You can use these inputs to simulate the following circuits:
%  1) AND(A,B)
%  2) AND4(A,B,C,D)
%  3) N-bit Carry_skip_adder(A,B,C), where the C is the carry_in
%  4) FA(A,B,C), where the C is the carry_in
%  5) HA(A,B)
%  6) Mux2to1(A,B,C), where the C is the selector
%  7) NOT(A)
%  8) OR(A,B)
%  9) N-bit RCA(A,B,C), where the C is the carry_in
%  10) XOR(A,B)


%%%%%%%%%%%%%%%%%%%%%%%% simulation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_bin = [0;1;1;0;1;1]; 
B_bin = [1;0;1;1;0;1]; 
C_bin = [1;1;1;0;0;0]; 
D_bin = [0];

Nbit = 2;   % it is used by RCA and CSA (parallelism). For the Carry-Skip Adder, the Nbit must be a multiple of 4, which is a constraint of the CSA model
opt_parameters = {'out_signal_plot','HA_without_regS','HA_without_regC','XOR_without_regS'}; % optional parameters, it can be empty
titleFontSize = 25;   % title FontSize of the plots
axisFontSize = 13;    % axes FontSize of the plots
labelFontSize = 20;   % labels FontSize of the plots
legendFontSize = 14;  % legend FontSize of the plots
line_width = 2;       % LineWidth of the lines  
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
circuit=0;  % initialization
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
analyzed_figures = 0;
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
        case 5
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            [HA_S,HA_C] = HA(A, B, model, opt_parameters{:})
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
        case 7
            A = DAC(A_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('"\n')
            NOT_out = NOT(A,model,opt_parameters{:})
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
        case 8
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            OR_out = OR(A,B,model,opt_parameters{:})
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
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
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
            
        case 10
            A = DAC(A_bin(ii,:),model);
            B = DAC(B_bin(ii,:),model);
            fprintf('Simulation %d:  A = "',ii)
            fprintf('%d',A_bin(ii,:))
            fprintf('", B = "')
            fprintf('%d',B_bin(ii,:))
            fprintf('"\n')
            XOR_out = XOR(A,B,model,opt_parameters{:})
            figures = findobj(0,'Type','Figure');
            New_figures = 0;
            for j=1:max(size(figures))-analyzed_figures % for each plot
                ax = figures(j).CurrentAxes;
                ax.Title.String = 'Simulation ' + string(ii);
                ax.Title.FontSize = titleFontSize;
                ax.Legend.FontSize = legendFontSize;
                ax.XAxis.FontSize = axisFontSize;
                ax.YAxis.FontSize = axisFontSize;
                ax.XLabel.FontSize = labelFontSize;
                ax.YLabel.FontSize = labelFontSize;
                New_figures = New_figures + 1;
            end
            analyzed_figures = analyzed_figures + New_figures;
    end
end
% LineWidth modification
lines = findobj('type','line');
for i=1:1:max(size(lines))
    lines(i).LineWidth = line_width;
end


end % if model~=4