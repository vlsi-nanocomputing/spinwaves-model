CMOS_45nm %DA CANCELLARE forse

% Source and drain diffusion length
lungh_diff=2.5.*fs_I; 		%um

% nMOS source and drain diffusion capacitance
Cbottom_n=lungh_diff.*Cj0n.*(1+Vdd./(2.*Pbn)).^(-Mjn);     		%pF/um

% nMOS source and drain sidewall capacitance
%Csidewall_n=1e6.*Cjswn.*(1+Vdd./(2.*Pbswn)).^(-Mswn); 			%pF./um
Csidewall_n=Cjswn.*(1+Vdd./(2.*Pbswn)).^(-Mswn);                %pF./um

% pMOS source and drain diffusion capacitance
Cbottom_p=lungh_diff.*Cj0p.*(1+Vdd./(2.*Pbp)).^(-Mjp); 			%pF./um

% pMOS source and drain sidewall capacitance
% Csidewall_p=1e6.*Cjswp.*(1+Vdd./(2.*Pbswp)).^(-Mswp);  			%pF./um
Csidewall_p=Cjswp.*(1+Vdd./(2.*Pbswp)).^(-Mswp);                %pF./um

%WD=fs_I;
Wn = Wgate_n;          %nm
Wp = beta*Wgate_n;     %nm
perim_N=2.*lungh_diff+Wn/1000;              %um
perim_P=2.*lungh_diff+2.*Wp/1000;           %um  % PERCHE' C'E' IL 2.*Wp?
% Junction capacitances for minimal nMOS and pMOS
Cjn =  (Cbottom_n.*Wn.*1e-3  + Csidewall_n.*perim_N);                   %pF
Cjp =  (Cbottom_p.*Wp.*1e-3  + Csidewall_p.*perim_P);                   %pF
%%Junction capacitances for two times minimal dimension nMOS and pMOS  
Cj2n = (Cbottom_n.*2*Wn.*1e-3  + Csidewall_n.*(2*lungh_diff + 2*Wn/1000)); %pF
Cj2p = (Cbottom_p.*2*Wp.*1e-3  + Csidewall_p.*(2*lungh_diff + 2*Wp)); %pF

%pA is the probability that input A is equal to 1
%pB is the probability that input b is equal to 1

%Wn = Wgate;          %nm
%Wp = 1.29*Wgate;     %nm

Igate_nmos = Igate_n*Wn.*1e-3;       %nA
Igate_pmos = Igate_p*Wp.*1e-3;       %nA

Ioff_nmos = Ioff_n*Wn.*1e-3;      %nA
Ioff_pmos = Ioff_p*Wp.*1e-3;      %nA

%P_static nand 2-input
P_stat_nand2 = Vdd*((1-pA)*(1-pB)*(Ioff_nmos + 2*Igate_pmos)+(1-pA)*pB*(2*Ioff_nmos + 2*Igate_nmos + Igate_pmos)+(pA*(1-pB)*(2*Ioff_nmos + Igate_pmos))+(pA*pB*(2*Ioff_pmos + 2*2*Igate_nmos)))*1e-9	;%W
Ioff_nd2 = Ioff_nmos*((1-pA)*(1-pB)+2*pB*(1-pA)+2*pA*(1-pB))+ 2*pA*pB*Ioff_pmos;        %%nA
Igate_nd2 = Igate_pmos*(2*(1-pA)*(1-pB)+pB*(1-pA)+pA*(1-pB))+ Igate_nmos*(2*pB*(1-pA)+4*pA*pB);         %%nA
%diff = Power_static_nd2 - Vdd*(Ioff_nd2+Igate_nd2)*1e-9;

%%output switching activity
alpha_a = pA*(1 - pA);
alpha_b = pB*(1 - pB);
Pout = 1 - pA*pB;
alpha_out = Pout*(1 - Pout);
%     Px1 = (pA*(1-pB))/(1-(1-pA)*(1-pB)); ANALIZZARE BENE QUESTA FORMULA!!
Px1 = ((pA*(1-pB)) + (1-pA)*(1-pB))/(1-(1-pA)*(1-pB));
alpha_x1 = Px1*(1-Px1);

%%capacitances
Coverlap_n=Wn/1000*Cgd0n; %pF
Coverlap_p=Wp/1000*Cgd0p; %pF
%Cin for one input in NAND2
Cin = Cox*Leff/1000*2*Wn/1000 + Cox*Leff/1000*Wp/1000 + 2*Coverlap_n +2*Coverlap_p;	%pF
% C interemediate node between the two nMOS
Cx1 = Cj2n;				%pF 
% Output capacitance
Cout = Cj2n + 2*Cjp;		%pF

%%total switching capacitance
Cs_nd2 = Cin*(alpha_a + alpha_b)+(Cx1*alpha_x1)+(Cout*alpha_out);	%pF
%%P_dyn nand 2-input
P_dyn_nand2 = 1/2 * Cs_nd2*1e-12 * f * Vdd^2 ;%W

%% INVERTER POWER

% static power inverter
P_stat_INV = Vdd*((1-pA)*(Ioff_nmos + Igate_pmos)+pA*(Ioff_pmos + Igate_nmos))*1e-9	;%W
Ioff_INV = Ioff_nmos*(1-pA)+ pA*Ioff_pmos;        %%nA
Igate_INV = Igate_pmos*(1-pA)+Igate_nmos*(1-pA);         %%nA
%diff = Power_static_INV - Vdd*(Ioff_INV+Igate_INV)*1e-9;

% dynamic Power inverter
Pout_INV=1-pA;
alpha_out_INV = Pout_INV*(1-Pout_INV);
Cin_INV = Cox*Leff/1000*2*Wn/1000 + Cox*Leff/1000*Wp/1000 + 2*Coverlap_n +2*Coverlap_p;	%pF
Cx1_INV= 0; %pF
Cout_INV= Cjn+Cjp;  %pF
Cs_INV = Cin_INV*(alpha_a)+(Cx1_INV*alpha_x1)+(Cout_INV*alpha_out_INV);
P_dyn_INV = 1/2 * Cs_INV*1e-12 * f* Vdd^2		;%W

%% Delay: NAND and INVERTER

%%delay nand 2-input
% % Cl = Cin;		%pF
% % Cnd2 = Cout + Cl;		%pF
% % Cf01 = Cox*Leff*1e-3*Wn + Cox*Leff*1e-3*Wp + Coverlap_n + Coverlap_p;	%pF	 Cin inverter
% % Cinv = Cjn + Cjp + Cf01;			%pF
% % Cmos = Cox*Leff*1e-3*Wn;	%pF
% % 
% % t_mos = (Cmos*Vdd/(Ion_n*Wn))*1e-12/1e-9;	%s		MOS delay	
% % t_inv = t_mos*Cinv/Cmos; %s		inverter delay
% % t_nd2 = t_inv*Cnd2/Cinv; %s		nand 2-input delay
% % f_max = 1/t_nd2; 		%Hz		max frequency nand 2-input
% % 
% % tpd_nand2 = t_nd2;        %s
% % tpd_inv= t_inv;                 %s


%%delay nand 2-input Elmore model
% tpdr_nand2 = log(2)*R*(Cout + Cg_FO1 +Cx); %GENERAL FORMULA

%          R = 1/(mueff_n*1e8*Cox*1e-12*(Wn/Leff)*1e-3*(Vdd-Vth_n));	%Ohm
R = 1/(mueff_n*1e8*Cox*1e-12*(Wn/Leff)*(Vdd-Vth_n));	%Ohm
% tpdr1 = ((2*mueff_p/mueff_n + 2)+4)*R*Cout*1e-12;		%s
tpdr_nand2 = log(2)*(Cout + Cin + Cx1)*R*1e-12; %just one Cin because is FO1
% tpdf1 = (2*Cout*1e-12)*R/2+(((2*mueff_p/mueff_n) + 2 + 4*1)*Cout*1e-12)*(R/2 + R/2);	%s
tpdf_nand2 =  log(2)*(Cx1*1e-12)*R/2+(Cout+Cin)*1e-12*(R/2 + R/2);	%s
tpd_nand2 = (tpdr_nand2 + tpdf_nand2)/2;	%s


%%delay INV Elmore model
tpdr_INV = log(2)*(Cin_INV + Cout_INV)*R*1e-12; %s
tpdf_INV =  tpdr_INV;	%s
tpd_INV = (tpdr_INV + tpdf_INV)/2;	%s


%% Area
A_INV = (Lgate/1000 * Wn/1000)+(Lgate/1000 * Wp/1000);            %um^2
A_nand2 = 2*(Lgate/1000 * 2*Wn/1000)+2*(Lgate/1000 * Wp/1000);  %um^2

%% To have the same name for Multiplier

% nand2
Ps_NAND2 = P_stat_nand2;  %W
Pdyn_NAND2 = P_dyn_nand2/f; %W*s
tau_NAND2 = tpd_nand2*1e9;     %ns
area_NAND2 = A_nand2; % um^2

% Inverter
Ps_INV = P_stat_INV;  %W
Pdyn_INV = P_dyn_INV/f; %W*s
tau_INV= tpd_INV*1e9;    %ns
area_INV = A_INV;  % um^2 



