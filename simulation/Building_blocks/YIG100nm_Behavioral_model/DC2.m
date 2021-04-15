function [out_S,out_C] = DC2(in_signal,model,varargin)

% This function describes the behavior of the ideal DC2 (without damping).
% It receives a signal, and gives 2 output signals(out_S,out_C).
% The function has one constraint:
% *) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad], delay[ns]]

SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;        % length of the coupling region  [nm]
gap2=10;        % the gap between the second coupled waveguides [nm]
limitation=limitation2;
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
% default values
disp_curves_flag = 0; % =1 Plot: dispersion curve. Display: Lc and pow_par
out_signal_plot_flag = 0; % =1 Plot: output signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= nargin-2   % -2 because the first 2 parameters are the required ones
    switch string(varargin{ii})
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
        case 'Lw'
            L2 = varargin{ii+1};
            ii = ii+1;
        case 'gap'
            gap2 = varargin{ii+1};
            ii = ii+1;
        case 'limitation'
            limitation = varargin{ii+1};
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


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%

d=w+gap2;       % [nm]
DC2_design = [h, w, d, B];
[wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC2 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC2_ff1=wm1./(2*pi);
DC2_ff2=wm2./(2*pi);

% nonlinear shift due to the input power
DC2_akx = in_signal(1);
DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(DC2_akx).^2;
DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(DC2_akx).^2;


DC2_ks = interp1(real(DC2_ff1_s),k1,in_signal(2));
DC2_kas = interp1(real(DC2_ff2_s),k1,in_signal(2));
DC2_Lc = pi/abs(DC2_ks-DC2_kas);  % [nm]
DC2_pow_par = cos(pi*L2/(2*DC2_Lc))^2; % [%]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep the same frequency and phase
out_S(2) = in_signal(2);
out_C(2) = in_signal(2);
out_S(3) = in_signal(3);
out_C(3) = in_signal(3);

out_S(1) = in_signal(1) * sqrt(DC2_pow_par);
out_C(1) = in_signal(1) * sqrt(1-DC2_pow_par);
out_C(1) = out_C(1)/sqrt(2); % damping 

% propagation delay
L_sing = 2*(5*h)/sin(0.3491);  % the length of the zone outside the coupled region
out_S(4) = in_signal(4) + DC_delay_calculation([L2, L_sing],model);
out_C(4) = out_S(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operations %%%%%%%%%%%%%%%%%%%%%%%

if disp_curves_flag == 1
    d = w+gap2; 
    DC2_design = [h, w, d, B];
    [wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
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
    
    DC2_ks = interp1(abs(DC2_ff1),k1,SW_frequency);
    DC2_kas = interp1(abs(DC2_ff2),k1,SW_frequency);
    DC2_Lc = pi/abs(DC2_ks-DC2_kas);  % [nm]
    DC2_pow_par = cos(pi*L2/(2*DC2_Lc))^2;
    fprintf('\n DC2: the Lc of the plot (dispersion curves) is: %dnm \n',DC2_Lc)
    fprintf('\n DC2: with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC2_Lc,DC2_pow_par*100)
end
    

if out_signal_plot_flag == 1
    signal_plotting([out_S;out_C],model,'DC2 out S','DC2 out C')
    fprintf('\n DC2: out_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_S(4),out_S(1),out_S(2),out_S(3),normalization(out_S(1),model))
    fprintf('\n DC2: out_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_C(4),out_C(1),out_C(2),out_C(3),normalization(out_C(1),model))
end


end
