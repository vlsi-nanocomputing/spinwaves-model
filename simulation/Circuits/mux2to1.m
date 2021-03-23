function [mux2to1_out] = mux2to1(in_0,in_1,sel,model,varargin)
% the multiplexer is implemented using the boolean function: 
% f = XS'+YS, with X and Y the inputs, and S the sel 


%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
if nargin == 5 % out_signal_plot
    if string(varargin{1}) == 'out_signal_plot'
        out_signal_plot_flag = 1;
    else
        error('Unsupported parameter: %s', string(varargin(1)))
    end
elseif nargin > 5
    error('Too many input arguments.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mux2to1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sel1,sel2] = duplicator(sel,model);
out1 = AND(in_0,NOT(sel1,model),model);
out2 = AND(in_1,sel2,model);
mux2to1_out = OR(out1,out2,model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(mux2to1_out,model,'mux2to1 out');
    fprintf('\n mux2to1: mux_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',mux2to1_out(4),mux2to1_out(1),mux2to1_out(2),mux2to1_out(3),normalization(mux2to1_out(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

