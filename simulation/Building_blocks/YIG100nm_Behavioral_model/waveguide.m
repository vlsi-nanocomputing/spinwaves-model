function [out_signal] = waveguide(in_signal,Lw,model,varargin)


SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%
h=30;           % thickness  [nm]
w=100;          % width  [nm]
B=0;            % external field [mT]
limitation = limitation_sing_waveg; % for single waveguide
akx = in_signal(1);
SW_frequency = in_signal(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
out_signal_plot_flag = 0;% =1 to plot and to display the output signal
disp_curves_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= nargin-3   % -3 because the first 3 parameters are the required ones
    switch string(varargin{ii})
        case 'dispersion_curves'
            disp_curves_flag = 1;
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
        case 'thickness' % we change the default DC parameter
            h = varargin{ii+1};
            ii = ii+1;
        case 'width'
            w = varargin{ii+1};
            ii = ii+1;
        case 'limitation'
            limitation = varargin{ii+1};
            ii = ii+1;
%         case 'group_velocity'
%             vgr_sing = varargin{ii+1};
%             ii = ii+1;
        case 'dx'
            dx = varargin{ii+1};
            ii = ii+1;
        case 'external_field'
            B = varargin{ii+1};
            ii = ii+1;
        otherwise
            error('Unsupported parameter: %s', string(varargin(ii)))
    end
    ii = ii + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
waveguide_design = [h, w, B];
[wm0, Tkx] = waveguide_equations(dkx, kmax, limitation, waveguide_design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% waveguide behavior %%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
ff0=wm0./(2*pi);

ff0_s = ff0+Tkx.*abs(akx).^2;
k0 = interp1(abs(ff0_s),k1,SW_frequency);  % wavenumber of the SW
delta_phase = k0*Lw; % [rad], phase shift accumulated along the waveguide 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude
out_signal(1) = akx;

% frequency
out_signal(2) = SW_frequency;

% phase
out_signal(3) = mod(in_signal(3) + delta_phase, pi);

% propagation delay
out_signal(4) = in_signal(4) + Lw/vgr_sing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operation %%%%%%%%%%%%%%%%%%%%%%%%
if out_signal_plot_flag == 1
    signal_plotting(out_signal,model,'waveguide out');
    fprintf('\n waveguide: out_signal = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d rad, normalized power = %d%% \n',out_signal(4),out_signal(1),out_signal(2),out_signal(3),normalization(out_signal(1),model))
end

if disp_curves_flag == 1
    N=size(k1);
    N=N(2);
    figure
    hold on
    plot(k1,real(ff0))
    plot(k1,SW_frequency*ones(1,N))
    hold off
    grid on
    xlabel('Wavenumber k [rad/nm]','FontSize',20)
    legend('Waveguide dispersion curve [GHz]','SW frequency [GHz]')
    axis([0 0.025 1.8 2.4])
    title('Dispersion curve of the waveguide ','FontSize',20)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
