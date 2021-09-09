function [out_S,out_C] = DC2(in_signal,model_parameters,plot_info,varargin)

% This function describes the behavior of the DC2 (with damping).
% It receives a signal, and gives 2 output signals(out_S,out_C).
% The function has one constraint:
% *) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad], delay [ns]]


%SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;        % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
B=0;            % external field [mT]
gap_region1=113; % the max gap of the region1 for the region1 discretization
gap_region3=210; % the max gap of the region1 for the region3 discretization
limitation = model_parameters.limitation2; % for gap=10nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
% default values
Lc_avg_flag = 0; % =1 Display: delta_phase, Lc_avg and pow_par
disp_curves_flag = 0; % =1 Plot: dispersion curve. Display: Lc and pow_par
out_signal_plot_flag = 1; % =1 Plot: output signals
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
            L2 = varargin{ii+1};
            ii = ii+1;
        case 'gap'
            gap2 = varargin{ii+1};
            ii = ii+1;
        case 'limitation'
            limitation = varargin{ii+1};
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_region1 = (gap_region1 - gap2) / sin(20*2*pi/360);  % [nm], length of region1
L_region3 = (gap_region3 - gap2) / sin(20*2*pi/360);  % [nm], length of region3
k1=model_parameters.dkx:model_parameters.dkx:model_parameters.kmax;
delta_phase = 0;  % phase accumulation due to the coupling
DC2_akx = in_signal(1); % amplitude of the input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DC2 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %%%%%%%%%%%% region 1 %%%%%%%%%%%%
dl = model_parameters.dx/2;  % L_region1 discretization resolution
dgap= -dl*sin(20*2*pi/360);
%gap: from 213nm(100+10+103) to 110nm(100+10)
N_cycle = ceil(L_region1/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle
       dl = L_region1 - (N_cycle-1)*dl;
       d = w + gap2;
    else
       d = w + gap_region1 + i1*dgap;
    end
    DC2_akx = DC2_akx*exp(-dl/model_parameters.x_freepath);  % losses
    DC2_design = [h, w, d, B];
    [wm1, wm2, Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC2_design);
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,model_parameters.SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,model_parameters.SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl;  % coupling contribution
end


                     %%%%%%%%%%%% region 2 %%%%%%%%%%%%
d = w + gap2; % constant gap in the coupled region
DC2_design = [h, w, d, B];
[wm1, wm2, Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC2_design);
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);
dl = model_parameters.dx;  % L_region2 discretization resolution
N_cycle = ceil(L2/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle 
        dl = L2 - (N_cycle-1)*dl;
    end
    DC2_akx = DC2_akx*exp(-dl/model_parameters.x_freepath); % losses
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,model_parameters.SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,model_parameters.SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end



                     %%%%%%%%%%%% region 3 %%%%%%%%%%%%
dl = model_parameters.dx/2; % L_region3 discretization resolution
dgap= dl*sin(20*2*pi/360);
%gap: from 110nm(100+10) to 310nm(100+10+200)
N_cycle = ceil(L_region3/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle
        dl = L_region3 - (N_cycle-1)*dl;
        d = d + dl*sin(20*2*pi/360);
    else
        d = w + gap2 + i1*dgap;
    end
    DC2_akx = DC2_akx*exp(-dl/model_parameters.x_freepath);
    DC2_design = [h, w, d, B];
    [wm1, wm2, Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC2_design);
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,model_parameters.SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,model_parameters.SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
end
Lc_avg = pi*L2/delta_phase;  % average Lc
DC2_pow_par = cos(pi*L2/(2*Lc_avg))^2; % DC2 power partition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep the same frequency and phase
out_S(2) = in_signal(2);
out_C(2) = in_signal(2);
out_S(3) = in_signal(3);
out_C(3) = in_signal(3);

% amplitudes
out_S(1) = DC2_akx * sqrt(DC2_pow_par);
out_C(1) = DC2_akx * sqrt(1-DC2_pow_par);


% propagation delay
L_sing = 2*(5*h)/sin(0.3491);  % the length of the zone outside the coupled region (length_region1 + length_region3)
out_S(4) = in_signal(4) + DC_delay_calculation([L2, L_sing],model_parameters);
out_C(4) = out_S(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operations %%%%%%%%%%%%%%%%%%%%%%%
if Lc_avg_flag == 1
    fprintf('\n DC2: the accumulated phase between two modes along the directional coupler 2 is %d rad \n',delta_phase)
    fprintf('\n DC2: the average coupling length Lc_avg of the directional coupler 2 is %dnm, and the number of jumps N = L2/Lc = %d \n',Lc_avg,L2/Lc_avg)
    fprintf('\n DC2: with Lc_avg = %dnm, Pout1/(Pout1+Pout2) = %d%% \n',Lc_avg,DC2_pow_par*100)
end

if disp_curves_flag == 1
    d = w+gap2; 
    DC2_design = [h, w, d, B];
    [wm1, wm2, DC2_Tkx] = DC_equations(model_parameters.dkx, model_parameters.kmax, limitation, DC2_design);
    DC2_ff1=wm1./(2*pi);
    DC2_ff2=wm2./(2*pi);
    
    N=size(k1);
    N=N(2);
    figure
    hold on
    plot(k1,real(DC2_ff1))
    plot(k1,real(DC2_ff2))
    plot(k1,SW_frequency*ones(1,N))
    hold off
    grid on
    xlabel('Wavenumber k [rad/nm]','FontSize',20)
    legend('DC2 symmetric mode [GHz]','DC2 Antisymmetric mode [GHz]','SW frequency [GHz]')
    axis([0 0.025 1.8 2.4])
    title('Dispersion curves of the DC2','FontSize',20)
    
    DC2_ks = interp1(abs(DC2_ff1),k1,model_parameters.SW_frequency);
    DC2_kas = interp1(abs(DC2_ff2),k1,model_parameters.SW_frequency);
    DC2_Lc = pi/abs(DC2_ks-DC2_kas);  % [nm]
    DC2_pow_par = cos(pi*L2/(2*DC2_Lc))^2;
    fprintf('\n DC2: the Lc of the plot (dispersion curves) is: %dnm \n',DC2_Lc)
    fprintf('\n DC2: with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC2_Lc,DC2_pow_par*100)
end
    

if out_signal_plot_flag == 1
    signal_plotting([out_S;out_C],model_parameters,'DC2 out S','DC2 out C')
    fprintf('\n DC2: out_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_S(4),out_S(1),out_S(2),out_S(3),normalization(out_S(1),model_parameters))
    fprintf('\n DC2: out_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_C(4),out_C(1),out_C(2),out_C(3),normalization(out_C(1),model_parameters))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end