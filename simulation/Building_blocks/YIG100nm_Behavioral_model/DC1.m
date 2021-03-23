function [out,out_I] = DC1(in_A,in_B,model,varargin)

% This function describes the behavior of the ideal DC1 (without damping).
% It receives 2 signals (A and B), and gives 2 output signals(out,out_I).
% The function has some constraints:
% 1) the 2 input signals have the same frequency
% 2) phase(B) - phase(A) = pi/2
% 3) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad]]

SW_parameters % script
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L1=370;         % length of the coupling region  [nm]
gap1=50;        % the gap between the coupled waveguides  [nm]
limitation=limitation1;
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
% default values
disp_curves_flag = 0; % =1  Plot: dispersion curve. Display: Lc and pow_par
out_signal_plot_flag = 0; % =1 Plot: output signals
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
        case 'Lw'
            L1 = varargin{ii+1};
            ii = ii+1;
        case 'gap'
            gap1 = varargin{ii+1};
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

d=w+gap1;       % [nm]
DC1_design = [h, w, d, B];
[wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC1 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1=dkx:dkx:kmax;
DC1_ff1=wm1./(2*pi);
DC1_ff2=wm2./(2*pi);

% nonlinear shift due to the input power
% DC1_akx = in_A(1);
% DC1_ff1_s = DC1_ff1 + DC1_Tkx.*abs(DC1_akx).^2;
% DC1_ff2_s = DC1_ff2 + DC1_Tkx.*abs(DC1_akx).^2;

DC1_ks = interp1(real(DC1_ff1),k1,in_A(2));   % [rad/nm]
DC1_kas = interp1(real(DC1_ff2),k1,in_A(2));  % [rad/nm]
DC1_Lc = pi/abs(DC1_ks-DC1_kas);         % [nm]
DC1_pow_par = cos(pi*L1/(2*DC1_Lc))^2;  %  [%]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%

% the single signals (in terms of amplitude) at the 2 outputs of the DC 
% the number 1 and 2 indicate the first and second (out I) outputs,
% respectively
in_A1_akx = in_A(1)*sqrt(DC1_pow_par);
in_A2_axk = in_A(1)*sqrt(1-DC1_pow_par);
in_B1_akx = in_B(1)*sqrt(1-DC1_pow_par);
in_B2_axk = in_B(1)*sqrt(DC1_pow_par);


if in_A(2) ~= in_B(2)
    display('DC1 not expected case: the input signals do not have the same frequency')
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
elseif (in_B(3)-in_A(3)) ~= pi/2
    display('DC1 not expected case: the input B is not shifted by pi/2 with respect to the intput A')
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
else
    out(1) = in_A1_akx+in_B1_akx; % constructive interference
    out(2) = in_A(2);
    out(3) = in_A(3);
    
    % destructive interference
    if in_B2_axk>in_A2_axk
        out_I(1) = in_B2_axk-in_A2_axk; % close to 0 
        out_I(2) = in_A(2);
        out_I(3) = in_B(3);  % pi/2
    else
        out_I(1) = in_A2_axk(1)-in_B2_axk(1); % close to 0
        out_I(2) = in_A(2);
        out_I(3) = in_A(3)-pi/2; % -pi/2
    end
end

%%% propagation delay
L_sing = 2*(5*h)/sin(0.3491);  % the length of the zone outside the coupled region
t_in = max(in_A(4), in_B(4));
out(4) = t_in + DC_delay_calculation([L1, L_sing],model);
out_I(4) = out(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional operations %%%%%%%%%%%%%%%%%%%%%%%

if disp_curves_flag == 1
    d = w+gap1; 
    DC1_design = [h, w, d, B];
    [wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);
    DC1_ff1=wm1./(2*pi);
    DC1_ff2=wm2./(2*pi);
    
    N=size(k1);
    N=N(2);
    figure
    hold on
    plot(k1,real(DC1_ff1))
    plot(k1,real(DC1_ff2))
    plot(k1,SW_frequency*ones(1,N))
    hold off
    grid on
    xlabel('Wavenumber k [rad/nm]','FontSize',20)
    legend('DC1 symmetric mode [GHz]','DC1 Antisymmetric mode [GHz]','SW frequency [GHz]')
    axis([0 0.025 1.8 2.4])
    title('Dispersion curves of the DC1','FontSize',20)
    
    DC1_ks = interp1(abs(DC1_ff1),k1,SW_frequency);
    DC1_kas = interp1(abs(DC1_ff2),k1,SW_frequency);
    DC1_Lc = pi/abs(DC1_ks-DC1_kas);  % [nm]
    DC1_pow_par = cos(pi*L1/(2*DC1_Lc))^2;
    fprintf('\n DC1: the Lc of the plot (dispersion curves) is %dnm \n',DC1_Lc)
    fprintf('\n DC1: with Lc = %dnm, Pout1/(Pout1+Pout2) = %d%%\n',DC1_Lc,DC1_pow_par*100)
end
    

if out_signal_plot_flag == 1
    signal_plotting([out;out_I],model,'DC1 output 1','DC1 output idle')
    fprintf('\n DC1: out1 = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out(4),out(1),out(2),out(3),normalization(out(1),model))
    fprintf('\n DC1: out_idle = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_I(4),out_I(1),out_I(2),out_I(3),normalization(out_I(1),model))
end


end

