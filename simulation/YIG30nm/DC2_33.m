% clear all
close all
% clc

% 5000 -> 0.031563207160646
% 2400 -> 0.015147446257228
SW_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10;           % thinckness  [nm]
w=30;          % width  [nm]
L2=2460;         % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
DC2_akx=6*0.015147446257228  % input SW of the DC2. Microwave field = 2mT 
SW_frequency = 2.29; % GHz
x_freepath = 9.07e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
delta_phase = 0;


limitation=limitation2;

DC1_design = [h, w, d, B];
[wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC1_design);

ff1 = wm1./(2*pi);
ff2 = wm2./(2*pi);
N=size(k1);
N=N(2);
hold on
plot(k1,real(ff1))
plot(k1,real(ff2))
plot(k1,SW_frequency*ones(1,N))
hold off
ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
ff2_s = ff2+Tkx.*abs(DC2_akx).^2;

DC1_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
DC1_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
delta_k = abs(DC1_ks-DC1_kas);
Lc = pi/delta_k
L2/Lc
% % % % % % % % % % len=300;
% % % % % % % % % % dl = len/100;
% % % % % % % % % % dgap= 103/100;
% % % % % % % % % % i1=1;
% % % % % % % % % % for l=dl:dl:len
% % % % % % % % % %     d = w + 113 - i1*dgap;
% % % % % % % % % %     DC2_akx = DC2_akx*exp(-dl/x_freepath);
% % % % % % % % % %     DC2_design = [h, w, d, B];
% % % % % % % % % %     cd common
% % % % % % % % % %     [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
% % % % % % % % % %     cd ..
% % % % % % % % % %     
% % % % % % % % % %     ff1=wm1./(2*pi);
% % % % % % % % % %     ff2=wm2./(2*pi);
% % % % % % % % % %     
% % % % % % % % % %     ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
% % % % % % % % % %     DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
% % % % % % % % % %     delta_k = abs(DC2_ks-DC2_kas); % rad/nm
% % % % % % % % % %     delta_phase = delta_phase + delta_k*dl; 
% % % % % % % % % %     i1=i1+1;
% % % % % % % % % % end
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % % d = w + 10;
% % % % % % % % % % DC2_design = [h, w, d, B];
% % % % % % % % % % cd common
% % % % % % % % % % [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
% % % % % % % % % % cd ..
% % % % % % % % % % ff1=wm1./(2*pi);
% % % % % % % % % % ff2=wm2./(2*pi);
% % % % % % % % % % dx = L2/100;
% % % % % % % % % % for i1=dx:dx:L2
% % % % % % % % % %     DC2_akx = DC2_akx*exp(-dx/x_freepath);
% % % % % % % % % %     ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
% % % % % % % % % %     DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
% % % % % % % % % %     delta_k = abs(DC2_ks-DC2_kas); % rad/nm
% % % % % % % % % %     delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
% % % % % % % % % % end
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % % % which has a length = 585 nm. 
% % % % % % % % % % 
% % % % % % % % % % len=585;
% % % % % % % % % % % Let us consider a finite partitioning of the 'len' into 100 subintervals
% % % % % % % % % % dl = len/100;
% % % % % % % % % % dgap=200/100;
% % % % % % % % % % i1=1;
% % % % % % % % % % for l=dl:dl:len
% % % % % % % % % %     d = w + 10 + i1*dgap;
% % % % % % % % % %     DC2_akx = DC2_akx*exp(-dl/x_freepath);
% % % % % % % % % %     DC2_design = [h, w, d, B];
% % % % % % % % % %     cd common
% % % % % % % % % %     [wm1, wm2, Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
% % % % % % % % % %     cd ..
% % % % % % % % % %     
% % % % % % % % % %     ff1=wm1./(2*pi);
% % % % % % % % % %     ff2=wm2./(2*pi);
% % % % % % % % % %     
% % % % % % % % % %     ff1_s = ff1+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     ff2_s = ff2+Tkx.*abs(DC2_akx).^2;
% % % % % % % % % %     DC2_ks = interp1(abs(ff1_s),k1,SW_frequency);  % rad/nm
% % % % % % % % % %     DC2_kas = interp1(abs(ff2_s),k1,SW_frequency); % rad/nm
% % % % % % % % % %     delta_k = abs(DC2_ks-DC2_kas); % rad/nm
% % % % % % % % % %     delta_phase = delta_phase + delta_k*dl; 
% % % % % % % % % %     
% % % % % % % % % %     i1=i1+1;
% % % % % % % % % % end


% Lc_avg_2mT = pi*L2/delta_phase % average Lc
% % Lc_avg_4mT = pi*L2/delta_phase % average Lc
% pow_par_2mT = cos(pi*L2/(2*Lc_avg_2mT))^2; % [%]
% pow_par_4mT = cos(pi*L2/(2*Lc_avg_4mT))^2; % [%]