% clear all
close all
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;         % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
DC2_akx=2*0.072290169039954/sqrt(2);  % input SW of the DC2. Microwave field = 2mT 
SW_frequency = 2.282; % GHz
x_freepath = 8.58e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
delta_phase = 0;
cd common
    SW_parameters
cd ..

limitation=0.63;

len=300;
dl = len/100;
dgap= 103/100;
i1=1;
for l=dl:dl:len
    d = w + 113 - i1*dgap;
    DC2_akx = DC2_akx*exp(-dl/x_freepath);
    DC2_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
    i1=i1+1;
end


d = w + 10;
DC2_design = [h, w, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
cd ..
ff1=wm1./(2*pi);
ff2=wm2./(2*pi);
dx = L2/100;
for i1=dx:dx:L2
    DC2_akx = DC2_akx*exp(-dx/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
end



% which has a length = 585 nm. 

len=585;
% Let us consider a finite partitioning of the 'len' into 100 subintervals
dl = len/100;
dgap=200/100;
i1=1;
for l=dl:dl:len
    d = w + 10 + i1*dgap;
    DC2_akx = DC2_akx*exp(-dl/x_freepath);
    DC2_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
    DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC2_ks-DC2_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
    
    i1=i1+1;
end


% Lc_avg_2mT = pi*L2/delta_phase % average Lc
% % Lc_avg_4mT = pi*L2/delta_phase % average Lc
% pow_par_2mT = cos(pi*L2/(2*Lc_avg_2mT))^2; % [%]
% pow_par_4mT = cos(pi*L2/(2*Lc_avg_4mT))^2; % [%]