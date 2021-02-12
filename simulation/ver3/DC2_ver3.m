function [out_S,out_C] = DC2_ver3(in_signal)

% This function describes the behavior of the ideal DC2 (without damping).
% It receives a signal, and gives 2 output signals(out_S,out_C).
% The function has one constraint:
% *) the input and output variables are vectors, and they are composed in
%    the following way:
%    [amplitude(dimensionless), frequency [GHz], phase [rad]]
cd common
SW_parameters % script
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;        % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
L_region1=300;  % [nm]
L_region3=585;  % [nm]
gap_region1=113;
gap_region3=210;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
delta_phase = 0;
cd common
    SW_parameters
cd ..
DC2_akx = in_signal(1);

                     %%%%%%%%%%%% region 1 %%%%%%%%%%%%
dl = dx/2;
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
    DC2_akx = DC2_akx*exp(-dl/x_freepath);
    DC2_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation2, DC2_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
end


                     %%%%%%%%%%%% region 2 %%%%%%%%%%%%
d = w + gap2;
DC2_design = [h, w, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation2, DC2_design);
cd ..
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);
dl = dx;
N_cycle = ceil(L2/dl);
for i1=1:1:N_cycle
    if i1 == N_cycle 
        dl = L2 - (N_cycle-1)*dl;
    end
    DC2_akx = DC2_akx*exp(-dl/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
end



                     %%%%%%%%%%%% region 3 %%%%%%%%%%%%
dl = dx/2;
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
    DC2_akx = DC2_akx*exp(-dl/x_freepath);
    DC2_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation2, DC2_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
end
Lc_avg = pi*L2/delta_phase;
DC2_pow_par = cos(pi*L2/(2*Lc_avg))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%
% we keep the same frequency and phase
out_S(2) = in_signal(2);
out_C(2) = in_signal(2);
out_S(3) = in_signal(3);
out_C(3) = in_signal(3);

out_S(1) = DC2_akx * sqrt(DC2_pow_par);
out_C(1) = DC2_akx * sqrt(1-DC2_pow_par);


out_S(4) = in_signal(4) + tpd_DC2;
out_C(4) = in_signal(4) + tpd_DC2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end