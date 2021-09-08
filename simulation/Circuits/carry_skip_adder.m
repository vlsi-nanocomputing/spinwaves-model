function [X] = carry_skip_adder(A,B,carry,Nbit,model, plot_info, varargin)
% N-bit carry skip adder
% Please use a N as a multiple of 4
% NB: the function describes the behavioral simulation of the device 
% (carry skip adder), the following structure is different with respect to
% the VHDL code structure.


%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 1;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_plot = 'no_plot';

if plot_info == "plot_all"
    sub_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW_parameters % script
N_block = Nbit/4;

%%% We need to duplicate A and B, where the first ones are sent to the
% FAs, and second ones are used to calculate the propagates.
A1 = zeros(Nbit,N_inf);
A2 = zeros(Nbit,N_inf);
B1 = zeros(Nbit,N_inf);
B2 = zeros(Nbit,N_inf);
for i=1:Nbit
    [A1(i,:),A2(i,:)] = duplicator(A(i,:),model,sub_plot, varargin{:});
    [B1(i,:),B2(i,:)] = duplicator(B(i,:),model,sub_plot, varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%% propagate calculation %%%%%%%%%%%%%%%%%%%%%%%%%
p = zeros(Nbit,N_inf); % propagates of every (ai,bi)
for i=1:Nbit
    p(i,:) = XOR( A1(i,:), B1(i,:) ,model,sub_plot, varargin{:} );
end
P_block = zeros(N_block,N_inf); % propagate of every block
for i=1:N_block
    P_block(i,:) = AND4( p((i-1)*4+1,:), p((i-1)*4+2,:), p((i-1)*4+3,:), p((i-1)*4+4,:), model,sub_plot, varargin{:} );
end


%%%%%%%%%%%%%%%%%%%%%%%% 4-bit RCA of every block %%%%%%%%%%%%%%%%%%%%%
if rem(Nbit,4) == 0
    for i=N_block:-1:1 % for every block
        RCA_out = RCA_Nbit( A2((i-1)*4+1:(i-1)*4+4, :), B2((i-1)*4+1:(i-1)*4+4, :), carry, 4 ,model,sub_plot, varargin{:}); % 4 is the number of FA of every RCA group
        X( (i-1)*4+2:(i-1)*4+5, :) = RCA_out(2:end, :);
        carry = mux2to1(RCA_out(1,:), carry, P_block(i,:) ,model,sub_plot, varargin{:});
    end
    X(1,:) = carry;
else
    display('please check your Nbit of the CSA, it must be a multiple of 4!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    i2 = Nbit+1;
    for i1 = 1:Nbit+1
        opt_par{i1} = 'CSA\_out' + string(i2-1);
        fprintf('\n %d-bit CSA: CSA_out%d = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',Nbit,i2-1,X(i1,4),X(i1,1),X(i1,2),X(i1,3),normalization(X(i1,1),model))
        i2 = i2-1;
    end
    signal_plotting(X,model,opt_par{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

