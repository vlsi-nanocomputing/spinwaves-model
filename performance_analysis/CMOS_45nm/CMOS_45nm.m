%%%Parameters 45 nm
%%%Parameters file
%%2005 --- LOP tech from TamTam  %%Start analysys selezionare questa
%%tecnologia e simulare ciò che vuoi simulare
%clear all; 
%clc 

%Da dove derivano Wgate_n, mueff_n

Vdd = 0.9; %V
Vth_n   = 0.29; %V      - EXT     presa dai calcoli di Vth_long
Vth_p   = -Vth_n; %V      - EXT     
epsilon_0 = 8.854188e-6 ; %pF/um   
epsilon_SiO2 = 3.9;
tox   = 0.0014 ; %um 
Cox = epsilon_0*epsilon_SiO2/tox;   %%pF/um^2
Lgate   = 45.1; %nm
Wgate_n = 90; %nm for nMOS 
Xj      = 20; %nm      - EXT
Cj0n    = 2.7e-3; %pF/um^2 %N-MOS junction capacitance per unit area under zero-bias conditions 
Cj0p    = 3.3e-3; %pF/um^2 % P-MOS junction capacitance per unit area under zero-bias conditions 
Cjswn   = 9.2e-4; %pF/um %N-MOS sidewall junction capacitance
Cjswp   = 8e-4; %pF/um % P-MOS sidewall junction capacitance 
Cgd0n   = 1.35e-4; %pF/um %N-MOS overlap capacitance between gate and drain per unit transistor width
Cgd0p   = 1e-4; %pF/um % P-MOS overlap capacitance between gate and drain per unit transistor width
Ion_n   = 543.1; %uA/um          presi da TamTam
Ion_p   = 292.89; %nA/um         presi da TamTam
Ioff_n  = 3.12; %nA/um       presi da TamTam
Ioff_p  = 3.12; %nA/um       presi da TamTam
Igate_n = 151; %nA/um       presi da TamTam
Igate_p = 47.8; %nA/um       presi da TamTam
mueff_n = 264.3780; %cm^2/Vs      - EXT *Arbitrario    calcolati con script Matlab mobility.m
mueff_p = 59.8150; %cm^2/Vs      - EXT *Arbitrario  calcolati con script Matlab mobility.m
beta = mueff_n/mueff_p;
Gamma = 0.8;        %Scaling factor for lateral diffusion
% Parametri di cui non si è capito il significato
Mjn     = 0.38; % [%]
Mjp     = 0.45; % [%]
Mswn    = 0.22; % [%]
Mswp    = 0.265; % [%]
Pbn     = 0.85; %V
Pbp     = 0.87; %V
Pbswn   = 0.67; %V
Pbswp   = 0.76; %V
pA      = 0.5; %Probability input signal A is equal to 1
pB      = 0.5; %Probability input signal B is equal to 1
f       = 500e6; %Hz         %%il gruppo di tecnologia lavora a questa frequenza per la loro cella


Leff=Lgate-Gamma*Xj;        %nm
fs_I=Leff./1000;            %%[Leff _ um]



