function [S,C] = FA(A,B,carry_in,model_parameters,plot_info,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 1;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
HA_varargin = varargin;
HA_plot = 'no_plot';
if plot_info == "plot_all"
    HA_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_info == "plot_all"
    fprintf("HA1");
end
[S1,C1] = HA(A,B,model_parameters,HA_plot,HA_varargin{:});
if plot_info == "plot_all"
    fprintf("HA2");
end
[S,C2] = HA(S1,carry_in,model_parameters,HA_plot,HA_varargin{:});
if plot_info == "plot_all"
    fprintf("HA3");
end
[C,~] = HA(C1,C2,model_parameters,HA_plot,HA_varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting([S;C],model_parameters,'FA sum','FA carry');
    fprintf('\n FA: FA_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',S(4),S(1),S(2),S(3),normalization(S(1),model_parameters))
    fprintf('\n FA: FA_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',C(4),C(1),C(2),C(3),normalization(C(1),model_parameters))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

