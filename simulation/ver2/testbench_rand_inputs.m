% test_bench 
clear all
clc


Nbit = 32;
N_simulation = 1; % number of simulations
err_rep_file = 'RCA_32bit_50sim.txt';
%%%%%%%%%%%%%%%%%%%%%% decimal input generation %%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(N_simulation,Nbit);
B = zeros(N_simulation,Nbit);
C = zeros(N_simulation,1);

for i = 1:N_simulation
    A_dec = randi([0,2^Nbit-1],1,1);
    B_dec = randi([0,2^Nbit-1],1,1);

%%%%%%%%%%%%%%%%%%%%%% digital input conversion %%%%%%%%%%%%%%%%%%%%%%%%
    A(i,:) = dec_to_bin(A_dec,Nbit);      % std_logic_vector(N downto 0)
    B(i,:) = dec_to_bin(B_dec,Nbit);      % std_logic_vector(N downto 0) 
    C(i,:) = randi([0,1],1,1);            % std_logic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% digital to analog conversion %%%%%%%%%%%%%%%%%%%%%%%

    in_A = DAC_ver2(A(i,:));
    in_B = DAC_ver2(B(i,:));
    in_C = DAC_ver2(C(i,:));

    % in_D = DAC_ver2(D);
    % in_E = DAC_ver2(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data-Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     output = RCA_Nbit_ver2(in_A,in_B,in_C,Nbit);
%     output = carry_skip_adder_ver2(in_A,in_B,in_C,Nbit);
%  	output = CLA_4bit_ver2(in_A,in_B,in_C);

    % in order to evaluate the error of result bits, in this point we keep
    % the amplitudes, and at the "max error evaluation" section we will
    % normalize them
    output_sig(i,:)= output(:,1); % amplitudes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% reference solution calculation %%%%%%%%%%%%%%%%%%%
    cd ../ver1
    exact_output(i,:) = RCA_Nbit_ver1(A(i,:),B(i,:),C(i,:),Nbit);
    cd ../ver2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% solutions comparison %%%%%%%%%%%%%%%%%%%%%%%
    % analog to digital conversion of the simulation solution
    output_bin(i,:) = ADC_ver2(output);      

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


%%%%%%%%%%%%%%%%%%%%%%%  error evaluation and report  %%%%%%%%%%%%%%%%%%
normalized_output = normalization(output_sig);

% f = fopen(err_rep_file,'w');
% for i= 1:N_simulation
%     for j = 1:Nbit+1
%         t = normalized_output(i,j);
%         if t >= 100
%             fprintf(f,'%3.4f  ',t);
%         elseif t >= 10
%             fprintf(f,'%3.4f   ',t);
%         else
%             fprintf(f,'%3.4f    ',t);
%         end
%     end
%     fprintf(f,'\n\n');
% end
% fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
