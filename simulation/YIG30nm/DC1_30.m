clear all
% close all
% clc

cd common
SW_parameters
cd ..
%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;           % thinckness  [nm]
w=30;           % width  [nm]
L1=260;         % length of the coupling region  [nm]
L_DC1=L1;       % length of the DC1 [nm], we currently use the same value of L1
gap1=10;        % the gap between the coupled waveguides  [nm]
d=w+gap1;       % [nm]
B=0;            % external field [mT]

gap_region1=50;
gap_region3=50;
L_region1 = (gap_region1 - gap1) / (2*sin(20*2*pi/360));  % [nm]
L_region3 = (gap_region3 - gap1) / (2*sin(20*2*pi/360));  % [nm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ak_A = SW_amplitude
ak_B = 0
ak_B = SW_amplitude

%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
limitation = limitation1;
delta_phase = 0;
                     %%%%%%%%%%%% region 1 %%%%%%%%%%%%
dl = dx/3;
dgap = -dl*sin(20*2*pi/360)*2;
%gap: from 550nm(100+450) to 150nm(100+50)
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
%gap: from 150nm(100+50) to 550nm(450+100)
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

Lc_avg = pi*L1/delta_phase
DC1_pow_par = cos(pi*L1/(2*Lc_avg))^2


