function [OR_out] = OR(in_A,in_B,model,plot_info, varargin)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,C] = HA(in_A,in_B,model,HA_plot, HA_varargin{:});
[OR_out,~] = HA(S,C,model,HA_plot, HA_varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(OR_out,model,'OR out');
    fprintf('\n OR: OR_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',OR_out(4),OR_out(1),OR_out(2),OR_out(3),normalization(OR_out(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
