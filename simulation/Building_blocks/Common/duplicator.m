function [out_signal1,out_signal2] = duplicator(in_signal,model_parameters,plot_info,varargin)
% The function receives an input signal, and it gives two same signals at
% the output.(The 3 signals are same)

%SW_parameters % script

in_signal2 = phase_shifter(DAC(0,model_parameters),pi/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
%%%%%% default values
out_signal_plot_flag = 1;% =1 to plot the out_S and the out_C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DC1_varargin,DC2_varargin,regS_varargin, regC_varargin, DC_without_regS_flag,DC_without_regC_flag] = decodeDCParameters(varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HA outputs plotting
DC1_plot = 'no_plot';

if plot_info == "plot_all"
    DC1_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[out_signal1, out_signal2] = DC1(in_signal, in_signal2,model_parameters,DC1_plot,DC1_varargin{:});

out_signal1 = amplifier(out_signal1,model_parameters.duplicator_gain1);
out_signal2 = amplifier(out_signal2,model_parameters.duplicator_gain2);
out_signal2 = phase_shifter(out_signal2,pi/2);


if out_signal_plot_flag == 1
    signal_plotting([out_signal1;out_signal2],model_parameters,'Out signal 1','Out signal 2');
    fprintf('\n Duplicator: Out signal 1 = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_signal1(4),out_signal1(1),out_signal1(2),out_signal1(3),normalization(out_signal1(1),model_parameters))
    fprintf('\n Duplicator: Out signal 2 = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_signal2(4),out_signal2(1),out_signal2(2),out_signal2(3),normalization(out_signal2(1),model_parameters))
end

end

