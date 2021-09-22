function [out_signal] = regenerator_C(in_signal,model_parameters,plot_info,varargin)

% This function describes the behavior of the regenerator (DC+amplifier)
% fot the carry bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad], delay [ns]]

%SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_out = 2.3;%4;%3.6;     % output amplifier gain
h=30;               % thinckness  [nm]
w=100;              % width  [nm]
Lw=1516;%1950;            % coupled length of the DC
gap=10;             % the gap between the second coupled waveguides [nm]
B=0;                % external field [mT]
limitation = model_parameters.limitation2; 
DC4_akx = in_signal(1); % amplitude of the input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
% default values
Lc_avg_flag = 0; % =1 Display: delta_phase, Lc_avg and pow_par
disp_curves_flag = 0; % =1 Plot: dispersion curve. Display: Lc and pow_par
out_signal_plot_flag = 1; % =1 Plot: output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= nargin-3   % -2 because the first 2 parameters are the required ones
    switch string(varargin{ii})
        case 'Lc_avg'
            Lc_avg_flag = 1;
        case 'dispersion_curves'
            disp_curves_flag = 1;
        case 'thickness'
            h = varargin{ii+1};
            ii = ii+1;
        case 'width'
            w = varargin{ii+1};
            ii = ii+1;
        case 'Lw'
            Lw = varargin{ii+1};
            ii = ii+1;
        case 'gap'
            gap = varargin{ii+1};
            ii = ii+1;
        case 'limitation'
            limitation = varargin{ii+1};
            ii = ii+1;
        case 'gain_out'
            gain_out = varargin{ii+1};
            ii = ii+1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_info == "no_plot"
    out_signal_plot_flag = 0;
end
%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = w + gap;       % [nm]
DC_design = [h, w, d, B];
[wm1, wm2, Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC_design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% regenerator operation %%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=model_parameters.dkx:model_parameters.dkx:model_parameters.kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);


delta_phase = 0;
dl=model_parameters.dx/2; % discretization resolution
N_cycle = ceil(Lw/dl);
for i1=1:N_cycle
    if i1 == N_cycle
        dl = Lw - (N_cycle-1)*dl;
    end
    DC4_akx = DC4_akx*exp(-dl/model_parameters.x_freepath);
    ff1_s = ff1+Tkx.*abs(DC4_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC4_akx).^2;
    DC4_ks = interp1(abs(ff1_s),k1,model_parameters.SW_frequency);  % rad/nm
    DC4_kas = interp1(abs(ff2_s),k1,model_parameters.SW_frequency); % rad/nm
    delta_k = abs(DC4_ks-DC4_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end


Lc_avg = pi*Lw/delta_phase;  % [nm], average Lc
pow_par = cos(pi*Lw/(2*Lc_avg))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = in_signal; % to have the same frequency and phase

% power splitting and amplification
out_signal(1) = DC4_akx * sqrt(pow_par);
out_signal = amplifier(out_signal,gain_out);


% propagation delay
out_signal(4) = in_signal(4) + DC_delay_calculation(Lw,model_parameters);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operations %%%%%%%%%%%%%%%%%%%%%%%
if Lc_avg_flag == 1
    fprintf('\n RegC: the accumulated phase between two modes along the regenerator ''C'' is %d rad \n',delta_phase)
    fprintf('\n RegC: the average coupling length Lc_avg of the regenerator ''C'' is %dnm, and the number of jumps N = Lw/Lc = %d \n',Lc_avg,Lw/Lc_avg)
    fprintf('\n RegC: with Lc_avg = %dnm, Pout1/(Pout1+Pout2) = %d%% \n',Lc_avg,pow_par*100)
end

if disp_curves_flag == 1
    d = w+gap; 
    DC_design = [h, w, d, B];
    [wm1, wm2, DC_Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC_design);
    DC_ff1=wm1./(2*pi);
    DC_ff2=wm2./(2*pi);
    
    N=size(k1);
    N=N(2);
    figure
    hold on
    plot(k1,real(DC_ff1))
    plot(k1,real(DC_ff2))
    plot(k1,SW_frequency*ones(1,N))
    hold off
    grid on
    xlabel('Wavenumber k [rad/nm]','FontSize',20)
    legend('Regenarator ''C'' symmetric mode [GHz]','Regenarator ''C'' Antisymmetric mode [GHz]','SW frequency [GHz]')
    axis([0 0.025 1.8 2.4])
    title('Dispersion curves of the regenarator ''C'' ','FontSize',20)
    
    DC_ks = interp1(abs(DC_ff1),k1,model_parameters.SW_frequency);
    DC_kas = interp1(abs(DC_ff2),k1,model_parameters.SW_frequency);
    DC_Lc = pi/abs(DC_ks-DC_kas);  % [nm]
    DC_pow_par = cos(pi*Lw/(2*DC_Lc))^2;
    fprintf('\n RegC: the Lc of the plot (dispersion curves) is %dnm \n',DC_Lc)
    fprintf('\n RegC: with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC_Lc,DC_pow_par*100)
end
    

if out_signal_plot_flag == 1
    signal_plotting(out_signal,model_parameters,'Regenerator ''C'' output')
    fprintf('\n RegC: out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_signal(4),out_signal(1),out_signal(2),out_signal(3),normalization(out_signal(1),model_parameters))
end

end

