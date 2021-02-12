
clear all


%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L1=370;         % length of the coupling region  [nm]
L_DC1=L1;       % length of the DC1 [nm], we currently use the same value of L1
gap1=50;        % the gap between the coupled waveguides  [nm]
d=w+gap1;       % [nm]
B=0;            % external field [mT]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
for limitation=0.53:0.01:0.61
    
limitation
delta_phase = 0;
cd common
    SW_parameters
cd ..


DC1_akx = sqrt(2)*0.0781;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DC1 operation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% region 1 %%%%%%%%%%%%%%%%%%%%%%%%%
len=585/4;
dl = len/100;
dgap = -100/100;
i1=1;
for l=dl:dl:len
    d = w + 150 + i1*dgap;
    DC1_akx = DC1_akx*exp(-dl/x_freepath);
    DC1_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC1_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC1_akx).^2;
    DC1_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
    i1=i1+1;
end
% display('region 1')
% delta_phase
% display('end')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% region 2 %%%%%%%%%%%%%%%%%%%%%%%%%
dx = L1/100;
d=w+gap1; 
DC1_design = [h, w, d, B];
cd common
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);
cd ..

ff1=wm1./(2*pi);
ff2=wm2./(2*pi);
for i1=dx:dx:L1
    DC1_akx = DC1_akx*exp(-dx/x_freepath);
    ff1_s = ff1+Tkx.*abs(DC1_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC1_akx).^2;
    DC1_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
end
% 
% display('region 2')
% delta_phase
% display('end')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% region 3 %%%%%%%%%%%%%%%%%%%%%%%%%
len=585/4;
dl = len/100;
dgap = 100/100;
i1=1;
for l=dl:dl:len
    d = w + 50 + i1*dgap;
    DC1_akx = DC1_akx*exp(-dl/x_freepath);
    DC1_design = [h, w, d, B];
    cd common
    [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);
    cd ..
    
    ff1=wm1./(2*pi);
    ff2=wm2./(2*pi);
    
    ff1_s = ff1+Tkx.*abs(DC1_akx).^2;
    ff2_s = ff2+Tkx.*abs(DC1_akx).^2;
    DC1_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
    DC1_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
    delta_k = abs(DC1_ks-DC1_kas); % rad/nm
    delta_phase = delta_phase + delta_k*dl; 
    
    i1=i1+1;
end


% display('region 3')
% delta_phase
% display('end')


Lc_avg = pi*L1/delta_phase
pow_par = cos(pi*L1/(2*Lc_avg))^2


end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     ff1_s = ff1+Tkx.*abs(DC1_akx).^2;
%     ff2_s = ff2+Tkx.*abs(DC1_akx).^2;
%     DC1_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
%     DC1_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
%     delta_k = abs(DC1_ks-DC1_kas); % rad/nm
%     Lc = pi/delta_k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% outputs generation %%%%%%%%%%%%%%%%%%%%%%%%%%

