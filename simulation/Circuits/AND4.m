function [AND4_out] = AND4(A,B,C,D,model_parameters,plot_info,varargin)
% AND port with 4 inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 1;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AND_plot = 'no_plot';
if plot_info == "plot_all"
    AND_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AND4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AND_out1 = AND(A,B,model_parameters,AND_plot,varargin{:});
AND_out2 = AND(C,D,model_parameters,AND_plot,varargin{:});
AND4_out = AND(AND_out1, AND_out2,model_parameters,AND_plot,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(AND4_out,model_parameters,'AND4 out');
    fprintf('\n AND4: AND4_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',AND4_out(4),AND4_out(1),AND4_out(2),AND4_out(3),normalization(AND4_out(1),model_parameters))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

