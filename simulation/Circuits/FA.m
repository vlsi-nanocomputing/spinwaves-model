function [S,C] = FA(A,B,carry_in,model,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
% if nargin == 5 % out_signal_plot
%     if string(varargin{1}) == 'out_signal_plot'
%         out_signal_plot_flag = 1;
%     else
%         error('Unsupported parameter: %s', string(varargin(1)))
%     end
% elseif nargin > 5
%     error('Too many input arguments.')
% end


HA_varargin = {};
ii=1;
while ii <= nargin-4   % -4 because the first 4 parameters are the required ones
    switch string(varargin{ii})
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
            if max(size(HA_varargin)) == 0
                HA_varargin{1} = 'out_signal_plot';
            else
                HA_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'HA_without_regS'
            if max(size(HA_varargin)) == 0
                HA_varargin{1} = 'HA_without_regS';
            else
                HA_varargin{end+1} = 'HA_without_regS';
            end
                
        case 'HA_without_regC'
            if max(size(HA_varargin)) == 0
                HA_varargin{1} = 'HA_without_regC';
            else
                HA_varargin{end+1} = 'HA_without_regC';
            end
            
        otherwise
            error('Unsupported parameter: %s', string(varargin(1)))
    end
    ii = ii + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("HA1");
[S1,C1] = HA(A,B,model,HA_varargin{:});
fprintf("HA2");
[S,C2] = HA(S1,carry_in,model,HA_varargin{:});
fprintf("HA3");
[C,~] = HA(C1,C2,model,HA_varargin{:});
% if xor_without_regs_flag == 1
%     C = XOR(C2,C2,model, 'XOR_without_regS');
% else
%     C = XOR(C1,C2,model);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting([S;C],model,'FA sum','FA carry');
    fprintf('\n FA: FA_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',S(4),S(1),S(2),S(3),normalization(S(1),model))
    fprintf('\n FA: FA_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',C(4),C(1),C(2),C(3),normalization(C(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

