function [mux2to1_out] = mux2to1(in_0,in_1,sel,model_parameters,plot_info, varargin)
% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 


%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 1;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
mux_plot = 'no_plot';
if plot_info == "plot_all"
    mux_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mux2to1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sel1,sel2] = duplicator(sel,model_parameters,mux_plot,varargin{:});
sel1_neg = NOT(sel1,model_parameters,mux_plot,varargin{:});
out1 = AND(in_0,sel1_neg,model_parameters,mux_plot,varargin{:});
out2 = AND(in_1,sel2,model_parameters,mux_plot,varargin{:});
mux2to1_out = OR(out1,out2,model_parameters,mux_plot,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(mux2to1_out,model_parameters,'mux2to1 out');
    fprintf('\n mux2to1: mux_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',mux2to1_out(4),mux2to1_out(1),mux2to1_out(2),mux2to1_out(3),normalization(mux2to1_out(1),model_parameters))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

