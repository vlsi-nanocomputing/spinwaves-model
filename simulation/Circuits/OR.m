function [OR_out] = OR(in_A,in_B,model,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
HA_varargin = {};
ii=1;
while ii <= nargin-3   % -3 because the first 3 parameters are the required ones
    switch string(varargin{ii})
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
            
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,C] = HA(in_A,in_B,model,HA_varargin{:});
[S,C] = HA(S,C,model,HA_varargin{:});
OR_out = S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(OR_out,model,'OR out');
    fprintf('\n OR: OR_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',OR_out(4),OR_out(1),OR_out(2),OR_out(3),normalization(OR_out(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
