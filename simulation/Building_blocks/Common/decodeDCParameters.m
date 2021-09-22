function [DC1_varargin,DC2_varargin,regS_varargin, regC_varargin, DC_without_regS_flag,DC_without_regC_flag] = decodeDCParameters(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameter flags %%%%%%%%%%%%%%%%%%%%%
%%%%%% default values
DC_without_regS_flag = 0; % =1 to replace the regS by an amplifier
DC_without_regC_flag = 0; % =1 to replace the regC by an amplifier
DC1_varargin = {};
DC2_varargin = {};
regS_varargin = {};
regC_varargin = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% optional parameters reception %%%%%%%%%%%%%%%%%%%%%%%
ii=1;
while ii <= length(varargin)   % -3 because the first 3 parameters are the required ones
    switch string(varargin{ii})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the DC1
        case 'DC1_Lc_avg'
            DC1_varargin{end+1} = 'Lc_avg';
            
        case 'DC1_dispersion_curves'
            DC1_varargin{end+1} = 'dispersion_curves';

            
        case 'DC1_out_signal_plot'
            DC1_varargin{end+1} = 'out_signal_plot';
            
        case 'DC1_thickness'
            DC1_varargin{end+1} = 'thickness';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'DC1_width'
            DC1_varargin{end+1} = 'width';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

            
        case 'DC1_Lw'
            DC1_varargin{end+1} = 'Lw';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

            
        case 'DC1_gap'
            DC1_varargin{end+1} = 'gap';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC1_limitation'
            DC1_varargin{end+1} = 'limitation';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC1_dx'
            DC1_varargin{end+1} = 'dx';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'DC1_external_field'
            DC1_varargin{end+1} = 'external_field';
            DC1_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the DC2
            
            
        case 'DC2_Lc_avg'
            DC2_varargin{end+1} = 'Lc_avg';
            
        case 'DC2_dispersion_curves'
            DC2_varargin{end+1} = 'dispersion_curves';
            
        case 'DC2_out_signal_plot'
            DC2_varargin{end+1} = 'out_signal_plot';
            
        case 'DC2_thickness'
            DC2_varargin{end+1} = 'thickness';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'DC2_width'
            DC2_varargin{end+1} = 'width';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC2_Lw'
            DC2_varargin{end+1} = 'Lw';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC2_gap'
            DC2_varargin{end+1} = 'gap';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC2_limitation'
            DC2_varargin{end+1} = 'limitation';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC2_dx'
            DC2_varargin{end+1} = 'dx';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'DC2_external_field'
            DC2_varargin{end+1} = 'external_field';
            DC2_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the regC
            
            
        case 'regC_Lc_avg'
                regC_varargin{end+1} = 'Lc_avg';
            
        case 'regC_dispersion_curves'
                regC_varargin{end+1} = 'dispersion_curves';
            
        case 'regC_out_signal_plot'
                regC_varargin{end+1} = 'out_signal_plot';
            
        case 'regC_thickness'
            regC_varargin{end+1} = 'thickness';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'regC_width'
            regC_varargin{end+1} = 'width';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_Lw'
            regC_varargin{end+1} = 'Lw';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_gap'
            regC_varargin{end+1} = 'gap';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_limitation'
            regC_varargin{end+1} = 'limitation';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_gain_out'
            regC_varargin{end+1} = 'gain_out';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_dx'
            regC_varargin{end+1} = 'dx';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regC_external_field'
            regC_varargin{end+1} = 'external_field';
            regC_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% optional parameters of the regS



        case 'regS_Lc_avg'
            regS_varargin{end+1} = 'Lc_avg';
            
        case 'regS_dispersion_curves'
            regS_varargin{end+1} = 'dispersion_curves';
            
        case 'regS_out_signal_plot'
            regS_varargin{end+1} = 'out_signal_plot';
            
        case 'regS_thickness'
            regS_varargin{end+1} = 'thickness';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;

        case 'regS_width'
            regS_varargin{end+1} = 'width';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_Lw1'
            regS_varargin{end+1} = 'Lw1';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_Lw2'
            regS_varargin{end+1} = 'Lw2';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_gap'
            regS_varargin{end+1} = 'gap';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_limitation'
            regS_varargin{end+1} = 'limitation';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_gain_in'
            regS_varargin{end+1} = 'gain_in';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_gain_interm'
            regS_varargin{end+1} = 'gain_interm';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_gain_out'
            regS_varargin{end+1} = 'gain_out';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_dx'
            regS_varargin{end+1} = 'dx';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'regS_external_field'
            regS_varargin{end+1} = 'external_field';
            regS_varargin{end+1} = varargin{ii+1};
            ii = ii+1;
            
        case 'DC_without_regS'
            DC_without_regS_flag = 1;           
        case 'DC_without_regC'
            DC_without_regC_flag = 1;  
            
        otherwise
            error('Unsupported parameter: %s', string(varargin(ii)))
    end
    ii = ii + 1;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

