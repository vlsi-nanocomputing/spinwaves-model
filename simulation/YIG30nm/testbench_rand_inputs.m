% test_bench 
clear all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%% simulation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbit = 32;
N_simulation = 1; % number of simulations
% result_rep_file = 'CSA_48bit_20sim.txt';
% result_rep_file = 'RCA_32bit_50sim.txt';
result_rep_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% decimal input generation %%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N_simulation,Nbit);
B = zeros(N_simulation,Nbit);
C = zeros(N_simulation,1);

for i = 1:N_simulation
    A_dec = randi([0,2^Nbit-1],1,1);
    B_dec = randi([0,2^Nbit-1],1,1);

%%%%%%%%%%%%%%%%%%%%%% decimal-binary conversion %%%%%%%%%%%%%%%%%%%%%%%%
    A(i,:) = dec_to_bin(A_dec,Nbit);      % std_logic_vector(N-1 downto 0)
    B(i,:) = dec_to_bin(B_dec,Nbit);      % std_logic_vector(N-1 downto 0) 
    C(i,:) = randi([0,1],1,1);            % std_logic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% digital-analog conversion %%%%%%%%%%%%%%%%%%%%%%%

    in_A = DAC(A(i,:));
    in_B = DAC(B(i,:));
    in_C = DAC(C(i,:));

    % in_D = DAC_ver2(D);
    % in_E = DAC_ver2(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data-Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    output = RCA_Nbit(in_A,in_B,in_C,Nbit);
%      output = carry_skip_adder(in_A,in_B,in_C,Nbit);
%  	output = CLA_4bit_ver2(in_A,in_B,in_C);

    % in order to evaluate the error of result bits, in this point we keep
    % the amplitudes, and at the "max error evaluation" section we will
    % normalize them
    output_sig(i,:)= output(:,1); % amplitudes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% reference solution calculation %%%%%%%%%%%%%%%%%%%
    cd ../ver1
    exact_output(i,:) = RCA_Nbit_ver1(A(i,:),B(i,:),C(i,:),Nbit);
    cd ../ver3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs comparison %%%%%%%%%%%%%%%%%%%%%%%
    % analog to digital conversion of the simulation solution
    output_bin(i,:) = ADC_ver3(output);      

end


correct = 0;
if output_bin == exact_output  % comparison
    display('the simulation result is correct')
    correct = 1;
else
    display('the simulation result is not correct')
end                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% error search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if correct == 0
    m = 1;
    for n=1:N_simulation % search of the errors
        if exact_output(n,:) == output_bin(n,:) 
            % there is some problem with ~=
        else
            % we need to decompose n into (i,j)
            err_A(m,:) = A(n,:);
            err_B(m,:) = A(n,:);
            m = m+1;
        end
    end
    [err_A err_B];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%  accuracy evaluation and report  %%%%%%%%%%%%%%%%%
normalized_output = normalization(output_sig);

if result_rep_flag == 1
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% result plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
for i=1:N_simulation
    scatter([1:1:Nbit+1], normalized_output(i,end:-1:1), 'filled')
end
axis([1, Nbit+1, -10, 110])
xlabel('bit','FontSize',20)
ylabel('Normalized power (%)','FontSize',20)
plot([1:1:Nbit+1],99*ones(1,Nbit+1))
plot([1:1:Nbit+1],101*ones(1,Nbit+1))
plot([1:1:Nbit+1],1*ones(1,Nbit+1))
plot([1:1:Nbit+1],-1*ones(1,Nbit+1))
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc