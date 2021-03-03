function [out,out_I] = DC1(in_A,in_B)

% This function describes the behavior of the DC1 (with damping).
% It receives 2 signals (A and B), and gives 2 output signals(out,out_I).
% The function has some constraints:
% 1) the 2 input signals have the same frequency
% 2) phase(B) - phase(A) = pi/2
% 3) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad], delay [ns]]

cd common
    SW_parameters
cd ..
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;           % thinckness  [nm]
w=30;          % width  [nm]
L1=260;         % length of the coupling region  [nm]
L_DC1=L1;       % length of the DC1 [nm], we currently use the same value of L1
gap1=20;        % the gap between the coupled waveguides  [nm]
B=0;            % external field [mT]

gap_region1=50; % [nm]
gap_region3=50; % [nm]
L_region1 = (gap_region1 - gap1) / (2*sin(20*2*pi/360));  % [nm]
L_region3 = (gap_region3 - gap1) / (2*sin(20*2*pi/360));  % [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
ak_A = in_A(1);
ak_B = in_B(1);
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
delta_phase = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC1 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %%%%%%%%%%%% region 1 %%%%%%%%%%%%
dl = dx/3;
dgap = -dl*sin(20*2*pi/360)*2;
%gap: from 80nm(30+50) to 40nm(30+10)
N_cycle = ceil(L_region1/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle
       dl = L_region1 - (N_cycle-1)*dl;
       d = w + gap1;
    else
       d = w + gap_region1 + i1*dgap;
    end
    
    ak_A = ak_A*exp(-dl/x_freepath);
    ak_B = ak_B*exp(-dl/x_freepath);
    DC1_design = [h, w, d, B];
    cd common
    [wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation1, DC1_design);
    cd ..
    
    DC1_ff1=wm1./(2*pi);
    DC1_ff2=wm2./(2*pi);
    
    DC1_ff1_s = DC1_ff1 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ff2_s = DC1_ff2 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ks = interp1(abs(DC1_ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(DC1_ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end


            %%%%%%%%%%%%% region 2 %%%%%%%%%%%%%%  
d = w+gap1; 
DC1_design = [h, w, d, B];
cd common
[wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation1, DC1_design);
cd ..
DC1_ff1=wm1./(2*pi);
DC1_ff2=wm2./(2*pi);

dl = dx/2;
N_cycle = ceil(L1/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle 
        dl = L1 - (N_cycle-1)*dl;
    end
    ak_A = ak_A*exp(-dl/x_freepath);
    ak_B = ak_B*exp(-dl/x_freepath);
    DC1_ff1_s = DC1_ff1 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ff2_s = DC1_ff2 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ks = interp1(abs(DC1_ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(DC1_ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
end



             %%%%%%%%%%%%%%%%%%% region 3 %%%%%%%%%%%%%%%%%
dl = dx/3;
dgap = dl*sin(20*2*pi/360)*2;
%gap: from 80nm(30+50) to 40nm(30+10)
N_cycle = ceil(L_region3/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle
        dl = L_region3 - (N_cycle-1)*dl;
        d = d + dl*sin(20*2*pi/360)*2;
    else
        d = w + gap1 + i1*dgap;
    end
    ak_A = ak_A*exp(-dl/x_freepath);
    ak_B = ak_B*exp(-dl/x_freepath);
    DC1_design = [h, w, d, B];
    cd common
    [wm1, wm2, DC1_Tkx] = DC_equations(dkx, kmax, limitation1, DC1_design);
    cd ..
    
    DC1_ff1=wm1./(2*pi);
    DC1_ff2=wm2./(2*pi);
    
    DC1_ff1_s = DC1_ff1 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ff2_s = DC1_ff2 + DC1_Tkx .* (abs(ak_A).^2 + abs(ak_B).^2);
    DC1_ks = interp1(abs(DC1_ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(DC1_ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
end

Lc_avg = pi*L1/delta_phase;
DC1_pow_par = cos(pi*L1/(2*Lc_avg))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%

% the single signals (in terms of amplitude) at the 2 outputs of the DC 
% the number 1 and 2 indicate the first and second (out I) outputs,
% respectively
in_A1_akx = ak_A*sqrt(DC1_pow_par);
in_A2_axk = ak_A*sqrt(1-DC1_pow_par);
in_B1_akx = ak_B*sqrt(1-DC1_pow_par);
in_B2_axk = ak_B*sqrt(DC1_pow_par);


if in_A(2) ~= in_B(2)
    display('DC1 unpredicted case: the input signals have not the same frequency')
    out = [inf, inf, inf, inf];
    out_I = [inf, inf, inf, inf];
elseif (in_B(3)-in_A(3)) ~= pi/2
    display('DC1 unpredicted case: the input B is not shifted by pi/2 with respect to the intput A')
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

% propagation delay
cd common
SW_parameters % script
cd ..
t_in = max(in_A(4), in_B(4));
out(4) = t_in + tpd_DC1;
out_I(4) = t_in + tpd_DC1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end