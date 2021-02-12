clear all
close all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%% parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=30;           % thinckness  [nm]
w=100;          % width  [nm]
L2=3000;         % length of the coupling region  [nm]
gap2=10;        % the gap between the coupled waveguides  [nm]
d=w+gap2;       % [nm]
B=0;            % external field [mT]
% DC2_akx=0.0138e-3/sqrt(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



DC2_design = [h, w, d, B];
%%%%%%%%%%%%%%%%%%%%% equations implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%
dkx=1e-3;
kmax=0.025;
k1=dkx:dkx:kmax;
% 
% limitation=0.74;
% limitation=0.595;
cd common
SW_parameters
cd ..



%%%%%%%%%%%%%%%%%%%%% prova: domanda 5
% a1 = 0.0138e-3/sqrt(2);
% a2 = 2*0.0138e-3/sqrt(2);
% i1=1;
% for limitation = 0.4:0.004:0.8
%     cd common
%     [wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
%     cd ..
%     DC2_ff1=wm1./(2*pi);
%     DC2_ff2=wm2./(2*pi);
%     
%     DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(a1).^2;
%     DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(a1).^2;    
%     DC2_ks = interp1(real(DC2_ff1_s),k1,2.282);
%     DC2_kas = interp1(real(DC2_ff2_s),k1,2.282);
%     Lc1(i1) = pi/abs(DC2_ks-DC2_kas);  % [nm]
%     
%     DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(a2).^2;
%     DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(a2).^2;    
%     DC2_ks = interp1(real(DC2_ff1_s),k1,2.282);
%     DC2_kas = interp1(real(DC2_ff2_s),k1,2.282);
%     Lc2(i1) = pi/abs(DC2_ks-DC2_kas);  % [nm]
%     
%     
%     i1=i1+1;
% end
% 
% limitation = [0.4:0.004:0.8];
% figure
% hold on
% plot(limitation, 214*ones(1,101),'LineWidth',2)
% plot(limitation, 230*ones(1,101),'LineWidth',2)
% plot(limitation, Lc1,'LineWidth',2)
% plot(limitation, Lc2,'LineWidth',2)
% hold off
% legend('214','230', 'Lc1', 'Lc2')
% % % % % % % % % i1=1;
% % % % % % % % % len=585;
% % % % % % % % % dl = len/100;
% % % % % % % % % dgap=200/100;
% % % % % % % % % for l=dl:dl:len
% % % % % % % % %     DC2_design = [h, w, d+i1*dgap, B];
% % % % % % % % %     cd common
% % % % % % % % %     [wm1_t, wm2_t, DC2_Tkx_t] = DC_equations(dkx, kmax, limitation, DC2_design);
% % % % % % % % %     cd ..
% % % % % % % % %     wm1 = wm1 + wm1_t;
% % % % % % % % %     wm2 = wm2 + wm2_t;
% % % % % % % % %     DC2_Tkx = DC2_Tkx + DC2_Tkx_t;
% % % % % % % % %     
% % % % % % % % %     i1=i1+1;
% % % % % % % % % end





% x=linspace(0,L2); % 99 intervals  [um]
% k1=dkx:dkx:kmax;
% DC2_ff1=wm1./(2*pi);
% DC2_ff2=wm2./(2*pi);

% %%%%%%%
% 1) cambiare il valore di akx, considerando kx in tutte le formule.
% 2) calcolare N volte DC2_ff1 e sommare tutti
%     for i = 1:1:N
%         DC2_ff1(i)++
%     end
% 3) e poi fare lo shift con damping
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_phase = 0;
% x_freepath = 8.58e3; % [nm], calculated from Single_dispersioncurve script
% dx = L2/100;  % sub-interval range, [nm]
% for i1=dx:dx:L2  % 100 cycles/intervals
%     DC2_akx = DC2_akx*exp(-dx/x_freepath);
%     DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(DC2_akx).^2;
%     DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(DC2_akx).^2;
%     DC2_ks = interp1(abs(DC2_ff1_s),k1,SW_frequency);  % rad/nm
%     DC2_kas = interp1(abs(DC2_ff2_s),k1,SW_frequency); % rad/nm
%     delta_k = abs(DC2_ks-DC2_kas); % rad/nm
%     delta_phase = delta_phase + delta_k*dx; % [rad], phase shift accumulated until this sub-interval
% end
% delta_phase1 = delta_phase;

% i1=1;
% len=585;
% dl = len/100;
% dgap=100/100;
% for l=dl:dl:len
%     DC2_design = [h, w, d+i1*dgap, B];
%     cd common
%     [wm1, wm2, DC2_Tkx] = DC_equations(dkx, kmax, limitation, DC2_design);
%     cd ..
%     DC2_ff1=wm1./(2*pi);
%     DC2_ff2=wm2./(2*pi);
%     
%     DC2_akx = DC2_akx*exp(-dl/x_freepath);
%     DC2_ff1_s = DC2_ff1+DC2_Tkx.*abs(DC2_akx).^2;
%     DC2_ff2_s = DC2_ff2+DC2_Tkx.*abs(DC2_akx).^2;
%     DC2_ks = interp1(abs(DC2_ff1_s),k1,SW_frequency);  % rad/nm
%     DC2_kas = interp1(abs(DC2_ff2_s),k1,SW_frequency); % rad/nm
%     delta_k = abs(DC2_ks-DC2_kas); % rad/nm
%     delta_phase = delta_phase + delta_k*dl; % [rad], phase shift accumulated until this sub-interval
%     
% %     i1=i1+1;
% end


% pi*L2/delta_phase



toc