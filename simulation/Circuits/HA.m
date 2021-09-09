function [out_S,out_C] = HA(in_A,in_B,model_parameters,plot_info,varargin)

in_B = phase_shifter(in_B,pi/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
%%%%%% default values
out_signal_plot_flag = 1;% =1 to plot the out_S and the out_C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DC1_varargin,DC2_varargin,regS_varargin, regC_varargin, DC_without_regS_flag,DC_without_regC_flag] = decodeDCParameters(varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HA outputs plotting
DC1_plot = 'no_plot';
DC2_plot = 'no_plot';
regS_plot = 'no_plot';
regC_plot = 'no_plot';

if plot_info == "plot_all"
    DC1_plot = 'plot_all';
    DC2_plot = 'plot_all';
    regS_plot = 'plot_all';
    regC_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% DC1 instantiation
[DC1_out,DC1_out_I] = DC1(in_A,in_B,model_parameters,DC1_plot,DC1_varargin{:});

% DC2 instantiation
[out_S,out_C] = DC2(DC1_out,model_parameters,DC2_plot,DC2_varargin{:});

%SW_parameters
% regS instantiation
if DC_without_regS_flag == 0
    out_S = regenerator_S(out_S,model_parameters,regS_plot,regS_varargin{:});
else    
    out_S = amplifier(out_S,gain_S);
end

% regC instantiation
if DC_without_regC_flag == 0
    out_C = regenerator_C(out_C,model_parameters,regC_plot,regC_varargin{:});
else
    out_C = amplifier(out_C,gain_C);
end

if out_signal_plot_flag == 1
    signal_plotting([out_S;out_C],model_parameters,'HA out S','HA out C');
    fprintf('\n HA: out_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_S(4),out_S(1),out_S(2),out_S(3),normalization(out_S(1),model_parameters))
    fprintf('\n HA: out_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_C(4),out_C(1),out_C(2),out_C(3),normalization(out_C(1),model_parameters))
end


end
