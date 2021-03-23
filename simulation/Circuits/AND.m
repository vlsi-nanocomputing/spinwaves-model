function [AND_out] = AND(in_A,in_B,model,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
if nargin == 4 % out_signal_plot
    if string(varargin{1}) == 'out_signal_plot'
        out_signal_plot_flag = 1;
    else
        error('Unsupported parameter: %s', string(varargin(1)))
    end
elseif nargin > 4
    error('Too many input arguments.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_B = phase_shifter(in_B,pi/2);
[DC1_out,DC1_out_I] = DC1(in_A,in_B,model);
[out_S, AND_out] = DC2(DC1_out,model);
AND_out = regenerator_C(AND_out,model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(AND_out,model,'AND out');
    fprintf('\n AND: AND_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',AND_out(4),AND_out(1),AND_out(2),AND_out(3),normalization(AND_out(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
