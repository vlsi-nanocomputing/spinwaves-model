% test_bench 
clear all
clc


Nbit = 4;
%%%%%%%%%%%%%%%%%%%%%% decimal input generation %%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:2^Nbit-1
for j = 0:2^Nbit-1
    
    A_dec = i;
    B_dec = j;


%%%%%%%%%%%%%%%%%%%%%% digital input setting %%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = dec_to_bin(A_dec,Nbit);      % std_logic_vector(N downto 0)
    B = dec_to_bin(B_dec,Nbit);      % std_logic_vector(N downto 0) 
    C = 0;                           % std_logic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% digital to analog conversion %%%%%%%%%%%%%%%%%%%%%%%

    in_A = DAC_ver2(A);
    in_B = DAC_ver2(B);
    in_C = DAC_ver2(C);

    % in_D = DAC_ver2(D);
    % in_E = DAC_ver2(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data-Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    output = RCA_Nbit_ver2(in_A,in_B,in_C,Nbit);
    
    % in order to evaluate the error of result bits, in this point we keep
    % the amplitudes, and at the "max error evaluation" section we will
    % normalize them
    output_sig(4*i+j+1,:)= output(:,1); % amplitudes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% reference solution calculation %%%%%%%%%%%%%%%%%%%
    exact_output(4*i+j+1,:) = RCA_Nbit_ver1(A,B,C,Nbit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% solutions comparison %%%%%%%%%%%%%%%%%%%%%%%
    % analog to digital conversion of the simulation solution
    output_bin(4*i+j+1,:) = ADC_ver2(output);      

end
end


correct = 0;
if output_bin == exact_output  % comparison
    display("the simulation result is correct")
    correct = 1;
else
    display("the simulation result is not correct")
end                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% error report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if correct == 0
    m = 1;
    for n=1:2^(2*Nbit) % search of the errors
        if exact_output(n,:) == output_bin(n,:) 
            % there is some problem with ~=
        else
            % we need to decompose n into (i,j)
            err_A(m,:) = dec_to_bin(fix((n-1)/4),Nbit);
            err_B(m,:) = dec_to_bin(rem(n-1,4),Nbit);
            m = m+1;
        end
    end
    [err_A err_B]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%  error evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%
normalized_output = (output_sig./0.153).^2*100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
