function [XOR_out] = XOR(in_A,in_B,model,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signals
xor_without_regs_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
% if nargin == 4 % out_signal_plot
%     if string(varargin{1}) == 'out_signal_plot'
%         out_signal_plot_flag = 1;
%     else
%         error('Unsupported parameter: %s', string(varargin(1)))
%     end
% elseif nargin > 4
%     error('Too many input arguments.')
% end

HA_varargin = {};
ii=1;
while ii <= nargin-3   % -3 because the first 3 parameters are the required ones
    switch string(varargin{ii})
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
            
        case 'XOR_without_regS'
            xor_without_regs_flag = 1;
            
        otherwise
            error('Unsupported parameter: %s', string(varargin(1)))
    end
    ii = ii + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_B = phase_shifter(in_B,pi/2);
[DC1_out,DC1_out_I] = DC1(in_A,in_B,model);
[out_S,out_C] = DC2(DC1_out,model);
if xor_without_regs_flag == 1
    XOR_out = out_S;
else
    XOR_out = regenerator_S(out_S,model);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(XOR_out,model,'XOR out');
    fprintf('\n AND: XOR_out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',XOR_out(4),XOR_out(1),XOR_out(2),XOR_out(3),normalization(XOR_out(1),model))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

