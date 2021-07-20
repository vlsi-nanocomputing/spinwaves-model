function [out_S,out_C] = HA(in_A,in_B,model,plot_info,varargin)

in_B = phase_shifter(in_B,pi/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
%%%%%% default values
DC1_opt_param_flag = 0; % =1 there are some optional parameters for the DC1
DC2_opt_param_flag = 0; % =1 there are some optional parameters for the DC2
regS_opt_param_flag = 0; % =1 there are some optional parameters for the regS
regC_opt_param_flag = 0; % =1 there are some optional parameters for the regC
out_signal_plot_flag = 1;% =1 to plot the out_S and the out_C
HA_without_regS_flag = 0; % =1 to replace the regS by an amplifier
HA_without_regC_flag = 0; % =1 to replace the regC by an amplifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= nargin-4   % -3 because the first 3 parameters are the required ones
    switch string(varargin{ii})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the DC1
        case 'DC1_Lc_avg'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'Lc_avg';
            else
                DC1_varargin{end+1} = 'Lc_avg';
            end
            
        case 'DC1_dispersion_curves'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'dispersion_curves';
            else
                DC1_varargin{end+1} = 'dispersion_curves';
            end
            
        case 'DC1_out_signal_plot'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'out_signal_plot';
            else
                DC1_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'DC1_thickness'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'thickness';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'thickness';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'DC1_width'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'width';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'width';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC1_Lw'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'Lw';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'Lw';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC1_gap'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'gap';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'gap';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC1_limitation'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'limitation';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'limitation';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC1_dx'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'dx';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'dx';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'DC1_external_field'
            if DC1_opt_param_flag == 0
                DC1_opt_param_flag = 1;
                DC1_varargin{1} = 'external_field';
                DC1_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC1_varargin{end+1} = 'external_field';
                DC1_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the DC2
            
            
        case 'DC2_Lc_avg'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'Lc_avg';
            else
                DC2_varargin{end+1} = 'Lc_avg';
            end
            
        case 'DC2_dispersion_curves'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'dispersion_curves';
            else
                DC2_varargin{end+1} = 'dispersion_curves';
            end
            
        case 'DC2_out_signal_plot'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'out_signal_plot';
            else
                DC2_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'DC2_thickness'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'thickness';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'thickness';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'DC2_width'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'width';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'width';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC2_Lw'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'Lw';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'Lw';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC2_gap'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'gap';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'gap';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC2_limitation'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'limitation';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'limitation';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'DC2_dx'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'dx';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'dx';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'DC2_external_field'
            if DC2_opt_param_flag == 0
                DC2_opt_param_flag = 1;
                DC2_varargin{1} = 'external_field';
                DC2_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                DC2_varargin{end+1} = 'external_field';
                DC2_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the regC
            
            
        case 'regC_Lc_avg'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'Lc_avg';
            else
                regC_varargin{end+1} = 'Lc_avg';
            end
            
        case 'regC_dispersion_curves'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'dispersion_curves';
            else
                regC_varargin{end+1} = 'dispersion_curves';
            end
            
        case 'regC_out_signal_plot'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'out_signal_plot';
            else
                regC_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'regC_thickness'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'thickness';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'thickness';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'regC_width'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'width';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'width';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_Lw'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'Lw';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'Lw';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_gap'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'gap';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'gap';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_limitation'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'limitation';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'limitation';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_gain_out'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'gain_out';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'gain_out';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_dx'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'dx';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'dx';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regC_external_field'
            if regC_opt_param_flag == 0
                regC_opt_param_flag = 1;
                regC_varargin{1} = 'external_field';
                regC_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regC_varargin{end+1} = 'external_field';
                regC_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the regS



        case 'regS_Lc_avg'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'Lc_avg';
            else
                regS_varargin{end+1} = 'Lc_avg';
            end
            
        case 'regS_dispersion_curves'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'dispersion_curves';
            else
                regS_varargin{end+1} = 'dispersion_curves';
            end
            
        case 'regS_out_signal_plot'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'out_signal_plot';
            else
                regS_varargin{end+1} = 'out_signal_plot';
            end
            
        case 'regS_thickness'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'thickness';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'thickness';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end

        case 'regS_width'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'width';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'width';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_Lw1'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'Lw1';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'Lw1';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_Lw2'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'Lw2';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'Lw2';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_gap'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'gap';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'gap';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_limitation'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'limitation';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'limitation';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_gain_in'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'gain_in';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'gain_in';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_gain_interm'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'gain_interm';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'gain_interm';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_gain_out'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'gain_out';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'gain_out';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_dx'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'dx';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'dx';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'regS_external_field'
            if regS_opt_param_flag == 0
                regS_opt_param_flag = 1;
                regS_varargin{1} = 'external_field';
                regS_varargin{2} = varargin{ii+1};
                ii = ii+1;
            else
                regS_varargin{end+1} = 'external_field';
                regS_varargin{end+1} = varargin{ii+1};
                ii = ii+1;
            end
            
        case 'HA_without_regS'
            HA_without_regS_flag = 1;           
        case 'HA_without_regC'
            HA_without_regC_flag = 1;  
            
        otherwise
            error('Unsupported parameter: %s', string(varargin(ii)))
    end
    ii = ii + 1;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HA outputs plotting
DC1_plot = 'no_plot';
DC2_plot = 'no_plot';
regS_plot = 'no_plot';
regC_plot = 'no_plot';

if plot_info == "plot_all"
    DC1_plot = 'plot_all';
    DC2_plot = 'plot_all';
    regS_plot = 'plot_all';
    regC_plot = 'plot_all';
elseif plot_info == "no_plot"
    out_signal_plot_flag = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% DC1 instantiation
if DC1_opt_param_flag == 0
    [DC1_out,DC1_out_I] = DC1(in_A,in_B,model,DC1_plot);
else
    [DC1_out,DC1_out_I] = DC1(in_A,in_B,model,DC1_plot,DC1_varargin{:});
end

% DC2 instantiation
if DC2_opt_param_flag == 0
    [out_S,out_C] = DC2(DC1_out,model,DC2_plot);
else
    [out_S,out_C] = DC2(DC1_out,model,DC2_plot,DC2_varargin{:});
end

SW_parameters
% regS instantiation
if HA_without_regS_flag == 0
    if regS_opt_param_flag == 0
        out_S = regenerator_S(out_S,model,regS_plot);
    else
        out_S = regenerator_S(out_S,model,regS_plot,regS_varargin{:});
    end
else    
    out_S = amplifier(out_S,gain_S);
end

% regC instantiation
if HA_without_regC_flag == 0
    if regC_opt_param_flag == 0
        out_C = regenerator_C(out_C,model,regC_plot);
    else
        out_C = regenerator_C(out_C,model,regC_plot,regC_varargin{:});
    end
else
    out_C = amplifier(out_C,gain_C);
end

if out_signal_plot_flag == 1
    signal_plotting([out_S;out_C],model,'HA out S','HA out C');
    fprintf('\n HA: out_S = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_S(4),out_S(1),out_S(2),out_S(3),normalization(out_S(1),model))
    fprintf('\n HA: out_C = u(t-t0) a sin(2 \x03c0 f t + \x03c6), where t0 = %d ns, a = %d, f = %d GHz and \x03c6 = %d, normalized power = %d%% \n',out_C(4),out_C(1),out_C(2),out_C(3),normalization(out_C(1),model))
end


end
