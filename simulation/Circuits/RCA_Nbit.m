function [X] = RCA_Nbit(A,B,carry,Nbit,model, varargin)
% N-bit ripple carry adder


%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
% if nargin == 6 % out_signal_plot
%     if string(varargin{1}) == 'out_signal_plot'
%         out_signal_plot_flag = 1;
%     else
%         error('Unsupported parameter: %s', string(varargin(1)))
%     end
% elseif nargin > 6
%     error('Too many input arguments.')
% end

RCA_varargin = {};
ii=1;
while ii <= nargin-5   % -5 because the first 5 parameters are the required ones
        switch string(varargin{ii})
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
            if max(size(RCA_varargin)) == 0
                RCA_varargin{1} = 'out_signal_plot';
            else
                RCA_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'HA_without_regS'
            if max(size(RCA_varargin)) == 0
                RCA_varargin{1} = 'HA_without_regS';
            else
                RCA_varargin{end+1} = 'HA_without_regS';
            end
                
        case 'HA_without_regC'
            if max(size(RCA_varargin)) == 0
                RCA_varargin{1} = 'HA_without_regC';
            else
                RCA_varargin{end+1} = 'HA_without_regC';
            end
            
        otherwise
            error('Unsupported parameter: %s', string(varargin(1)))
    end
    ii = ii + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW_parameters % script
X = zeros(Nbit+1,N_inf);

for j=0:Nbit-1
    [X(Nbit-j+1,:),carry] = FA( A(Nbit-j,:), B(Nbit-j,:), carry, model, RCA_varargin{:});
end
X(1,:) = carry; % MSB or final carry_out of the RCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    i2 = Nbit+1;
    for i1 = 1:Nbit+1
        opt_par{i1} = 'RCA\_out' + string(i2-1);
        fprintf('\n %d-bit RCA: RCA_out%d = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',Nbit,i2-1,X(i1,4),X(i1,1),X(i1,2),X(i1,3),normalization(X(i1,1),model))
        i2 = i2-1;
    end
    signal_plotting(X,model,opt_par{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end