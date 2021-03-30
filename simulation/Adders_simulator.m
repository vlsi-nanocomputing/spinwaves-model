% Adders simulator
clear all % do not remove
tic

%%%%%%%%%%%%%%%%%%%%%%%% simulation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbit = 2;  % parallelism of the adder. Note that if you choose the CSA, the Nbit must be a mulple of 4, that is a constraint of the CSA model
N_simulation = 1; % number of simulations
result_rep_file = 'RCA_8bit_5sim.txt'; % result report file_name 
result_rep_flag = 0; % =1 to write the output powers (in terms of normalized power) in "the result_rep_file".txt
err_search_flag = 0; % =1 to search and to display the combinations that generated wrong outputs
plot_range_flag = 0; % =1 to set ranges for logic '1' and '0' in the final plot, using the following parameters:
logic1_sup = 101;
logic1_inf = 99;
logic0_sup = 1;
logic0_inf = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Choosing of model category %%%%%%%%%%%%%%%%%%%
model=0;
while model ~= [1,2,3,4]
    model = input('Choose one model category from the following list:\n 1) YIG (100 nm) Behavioral Model \n 2) YIG (100 nm) Physical Model \n 3) YIG (30 nm) Physical Model \n 4) QUIT \n');
end
if model~=4  % if model==4, the program terminates

switch model % to get the correct folder for building blocks
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Choosing of the simulation adder %%%%%%%%%%%%%%%%
adder=0;
while adder ~= [1,2]
    adder = input('\nChoose one simulation adder from the following list: \n  1) Ripple-Carry Adder \n  2) Carry-Skip Adder \n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% decimal input generation %%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N_simulation,Nbit);
B = zeros(N_simulation,Nbit);
C = zeros(N_simulation,1);

for i = 1:N_simulation 
    A_dec = randi([0,2^Nbit-1],1,1); % random inputs generation
    B_dec = randi([0,2^Nbit-1],1,1);

%%%%%%%%%%%%%%%%%%%%%% decimal-binary conversion %%%%%%%%%%%%%%%%%%%%%%%%
    
    A(i,:) = dec_to_bin(A_dec,Nbit);      % std_logic_vector(N-1 downto 0)
    B(i,:) = dec_to_bin(B_dec,Nbit);      % std_logic_vector(N-1 downto 0) 
    C(i,:) = randi([0,1],1,1);            % std_logic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% digital-analog conversion %%%%%%%%%%%%%%%%%%%%%%%
    in_A = DAC(A(i,:),model);
    in_B = DAC(B(i,:),model);
    in_C = DAC(C(i,:),model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% you can choose one adder from the following list:

    if adder == 1
        output = RCA_Nbit(in_A,in_B,in_C,Nbit,model);
    elseif adder == 2
        output = carry_skip_adder(in_A, in_B, in_C, Nbit, model);
    end
    
    % in order to evaluate the error of result bits, in this point we store
    % their amplitudes, and at the "accuracy evaluation and report" section we will
    % normalize them
    output_sig(i,:)= output(:,1)'; % amplitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% reference solution calculation %%%%%%%%%%%%%%%%%%%
    model_t = model;
    model = 0; % logical model
    addpath('Building_blocks/YIG100nm_Logical_model')% model category change
    exact_output(i,:) = RCA_Nbit(A(i,:)',B(i,:)',C(i,:),Nbit,model);
    addpath(model_path)
    model = model_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% outputs analog-digital conversion %%%%%%%%%%%%%%%%%
    % analog to digital conversion of the simulation solution
    output_bin(i,:) = ADC(output,model);      
end % for i = 1:N_simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs comparison %%%%%%%%%%%%%%%%%%%%
correct = 0;
if output_bin == exact_output  % comparison
    display('The simulation result is correct.')
    correct = 1;
else
    display('The simulation result is not correct.')
end                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% error search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if err_search_flag == 1 % you can set this flag at the beginning
    if correct == 0
        m = 1;
        for n=1:N_simulation % search of the errors
            if exact_output(n,:) ~= output_bin(n,:) 
                % we need to decompose n into (i,j)
                err_A(m,:) = A(n,:);
                err_B(m,:) = B(n,:); % A and B are std_logic_vectors
                err_C(m,:) = C(n,:);
                m = m+1;
            end
        end
        display('The following combinations (A and B) generated wrong outputs:')
        m = m-1;
        for i=1:1:m % error combination display
            fprintf('\nCombination %d: \n \t A = ',i)
            fprintf('%d',err_A(i,:))
            fprintf('\n \t B = ')
            fprintf('%d',err_B(i,:))
            fprintf('\n \t C = ')
            fprintf('%d',err_C(i,:))
        end
        fprintf('\n')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%  accuracy evaluation and storing  %%%%%%%%%%%%%%%%%
% normalization of the output bits of every simulation
normalized_output = normalization(output_sig,model); 

if result_rep_flag == 1 % report of the simulation result
    f = fopen(result_rep_file,'w');
    for i= 1:N_simulation
        for j = 1:Nbit+1
            t = normalized_output(i,j);
            if t >= 100
                fprintf(f,'%3.4f  ',t);
            elseif t >= 10
                fprintf(f,'%3.4f   ',t);
            else
                fprintf(f,'%3.4f    ',t);
            end
        end
        fprintf(f,'\n\n');
    end
    fclose(f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% accuracy plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
for i=1:N_simulation
    scatter([1:1:Nbit+1], normalized_output(i,end:-1:1), 'filled')
end
axis([1, Nbit+1, -10, 110])
xlabel('bit','FontSize',20)
ylabel('Normalized power (%)','FontSize',20)
if plot_range_flag == 1
    plot([1:1:Nbit+1],logic1_sup*ones(1,Nbit+1))
    plot([1:1:Nbit+1],logic1_inf*ones(1,Nbit+1))
    plot([1:1:Nbit+1],logic0_sup*ones(1,Nbit+1))
    plot([1:1:Nbit+1],logic0_inf*ones(1,Nbit+1))
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
toc