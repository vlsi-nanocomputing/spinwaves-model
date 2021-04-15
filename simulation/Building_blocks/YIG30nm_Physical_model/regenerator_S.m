function [out_signal] = regenerator_S(in_signal,model,varargin)

% This function describes the behavior of the regenerator (2DC+2amplifier)
% for the sum bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad], delay [ns]]

SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_in = 8;        % input amplifier gain
gain_interm = 2;    % intermediate amplifier gain
h=10;               % thinckness  [nm]
w=30;               % width  [nm]
L3=900;             % length of the first DC coupling region  [nm]
L4=978;             % length of the second DC coupling region  [nm]
gap=10;             % the gap between the second coupled waveguides [nm]
B=0;                % external field [mT]
limitation = limitation2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
% default values
Lc_avg_flag = 0; % =1 Display: delta_phase, Lc_avg and pow_par
disp_curves_flag = 0; % =1 Plot: dispersion curve. Display: Lc and pow_par
out_signal_plot_flag = 0; % =1 Plot: output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= nargin-2   % -2 because the first 2 parameters are the required ones
    switch string(varargin{ii})
        case 'Lc_avg'
            Lc_avg_flag = 1;
        case 'dispersion_curves'
            disp_curves_flag = 1;
        case 'out_signal_plot'
            out_signal_plot_flag = 1;
        case 'thickness'
            h = varargin{ii+1};
            ii = ii+1;
        case 'width'
            w = varargin{ii+1};
            ii = ii+1;
        case 'Lw1'
            L3 = varargin{ii+1};
            ii = ii+1;
        case 'Lw2'
            L4 = varargin{ii+1};
            ii = ii+1;
        case 'gap'
            gap = varargin{ii+1};
            ii = ii+1;
        case 'limitation'
            limitation = varargin{ii+1};
            ii = ii+1;
        case 'gain_in'
            gain_in = varargin{ii+1};
            ii = ii+1;
        case 'gain_interm'
            gain_interm = varargin{ii+1};
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

%%%%%%%%%%%%%%%%%%%%% input signal amplification %%%%%%%%%%%%%%%%%%%%%%%%%
DC3_akx = amplifier(in_signal,gain_in);
DC3_akx = DC3_akx(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
d=w+gap;       % [nm]
DC_design = [h, w, d, B];
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC_design);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% first DC operation %%%%%%%%%%%%%%%%%%%%%%%%%%%

k1=dkx:dkx:kmax;
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);

delta_phase = 0;
dl=dx/2;
N_cycle = ceil(L3/dl);
for i1=1:N_cycle
    if i1 == N_cycle
        dl = L3 - (N_cycle-1)*dl;
    end
    DC3_akx = DC3_akx*exp(-dl/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC3_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC3_akx).^2;
    DC3_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC3_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC3_ks-DC3_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end


Lc_avg = pi*L3/delta_phase;  % [nm]
pow_par = cos(pi*L3/(2*Lc_avg))^2;
DC4_akx = DC3_akx * sqrt(pow_par);
DC4_akx = [DC4_akx, in_signal(2:4)]; % I need the complete SW vector for the amplifier function
DC4_akx = amplifier(DC4_akx,gain_interm);
DC4_akx = DC4_akx(1);
%%%%% these parameters can be used in the "optional operations" section %%%
delta_phase1 = delta_phase;
Lc_avg1 = Lc_avg;
pow_par1 = pow_par;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% second DC operation %%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_phase = 0;
dl=dx/2;
N_cycle = ceil(L4/dl);
for i1=1:N_cycle
    if i1 == N_cycle
        dl = L4 - (N_cycle-1)*dl;
    end
    DC4_akx = DC4_akx*exp(-dl/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC4_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC4_akx).^2;
    DC4_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC4_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC4_ks-DC4_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end

Lc_avg = pi*L4/delta_phase;  % [nm], average Lc
pow_par = cos(pi*L4/(2*Lc_avg))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
out_signal = in_signal; % to have the same frequency and phase


% power splitting and amplification
out_signal(1) = DC4_akx * sqrt(1-pow_par);

% propagation delay
out_signal(4) = in_signal(4) + DC_delay_calculation(L3+L4,model); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operations %%%%%%%%%%%%%%%%%%%%%%%
if Lc_avg_flag == 1
    fprintf('\n RegS: the accumulated phases between two modes along the two DCs are: %d rad and %d rad \n',delta_phase1,delta_phase)
    fprintf('\n RegS: the average coupling length Lc_avg of the first DC is %dnm, and the number of jumps N = Lw/Lc = %d \n',Lc_avg1,L3/Lc_avg1)
    fprintf('\n RegS: the average coupling length Lc_avg of the second DC is %dnm, and the number of jumps N = Lw/Lc = %d \n',Lc_avg,L4/Lc_avg)
    fprintf('\n RegS: for the first DC output, with average Lc_avg = %dnm, Pout1/(Pout1+Pout2) = %d%% \n',Lc_avg1,pow_par1*100)
    fprintf('\n RegS: for the second DC output, with average Lc_avg = %dnm, Pout1/(Pout1+Pout2) = %d%% \n',Lc_avg,pow_par*100)
end

if disp_curves_flag == 1
    d = w+gap; 
    DC_design = [h, w, d, B];
    [wm1, wm2, DC_Tkx] = DC_equations(dkx, kmax, limitation, DC_design);
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
    legend('Regenarator ''S'' symmetric mode [GHz]','Regenarator ''S'' Antisymmetric mode [GHz]','SW frequency [GHz]')
    axis([0 0.025 1.8 2.4])
    title('Dispersion curves of the regenarator ''S'' (for both the DCs)','FontSize',15)
    
    DC_ks = interp1(abs(DC_ff1),k1,SW_frequency);
    DC_kas = interp1(abs(DC_ff2),k1,SW_frequency);
    DC_Lc = pi/abs(DC_ks-DC_kas);  % [nm]
    fprintf('\n RegS: the Lc of the plot (dispersion curves) is %dnm \n',DC_Lc)
    DC_pow_par = cos(pi*L3/(2*DC_Lc))^2;
    fprintf('\n RegS: for the first DC, with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC_Lc,DC_pow_par*100)
    DC_pow_par = cos(pi*L4/(2*DC_Lc))^2;
    fprintf('\n RegS: for the second DC, with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC_Lc,DC_pow_par*100)
end
    

if out_signal_plot_flag == 1
    signal_plotting(out_signal,model,'Regenerator ''S'' output')
    fprintf('\n RegS: out = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_signal(4),out_signal(1),out_signal(2),out_signal(3),normalization(out_signal(1),model))
end


end

