function [tpd_DC] = DC_delay_calculation(lengths,model)


switch model % choosing of vgr according to the technology
    case 1  % YIG100nm behavioral model
        vgr_sing = 195.3;    % group velocity of single waveguide [m/s]
        vgr_coup = 25;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
    case 2 % YIG100nm physical model
        vgr_sing = 195.3;    % group velocity of single waveguide [m/s]
        vgr_coup = 25;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
    case 3 % YIG30nm physical model
        vgr_sing = 205.7;    % group velocity of single waveguide [m/s]
        vgr_coup = 137;       % group velocity in the coupled region, that is the vgr of antisymmetric mode [m/s]
end

% type of DC
DC_type = size(lengths);
DC_type = DC_type(2); % =2 it's DC1 or DC2, =1 it's regS or regC
% lengths(1) is the length of coupling region, lengths(2) is the other
% regions one

if DC_type == 1
    tpd_DC = lengths/vgr_coup;      % [ns]
else
    tpd_DC = lengths(1)/vgr_coup + lengths(2)/vgr_sing;        % [ns] 
end

end

