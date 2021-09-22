function [AND_out] = AND(in_A,in_B,model_parameters,plot_info,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 1;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
[DC1_varargin,DC2_varargin,regS_varargin, regC_varargin, DC_without_regS_flag,DC_without_regC_flag] = decodeDCParameters(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DC1_plot = 'no_plot';
DC2_plot = 'no_plot';
regC_plot = 'no_plot';

if plot_info == "plot_all"
    DC1_plot = 'plot_all';
    DC2_plot = 'plot_all';
    regC_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_B = phase_shifter(in_B,pi/2);
[DC1_out,DC1_out_I] = DC1(in_A,in_B,model_parameters,DC1_plot,DC1_varargin{:});
[out_S, AND_out] = DC2(DC1_out,model_parameters,DC2_plot,DC2_varargin{:});
% regC instantiation
if DC_without_regC_flag == 0
    AND_out = regenerator_C(AND_out,model_parameters,regC_plot,regC_varargin{:});
else
    AND_out = amplifier(out_C,gain_C);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(AND_out,model_parameters,'AND out');
    fprintf('\n AND: AND_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',AND_out(4),AND_out(1),AND_out(2),AND_out(3),normalization(AND_out(1),model_parameters))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
