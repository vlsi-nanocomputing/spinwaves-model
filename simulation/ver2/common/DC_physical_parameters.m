% This script contains all common physical parameters of directional coupler 


%%%%%%%%%%%%%%%%%% physical parameters (constants) %%%%%%%%%%%%%%%%%%%%
Ms=1.4e5;       % Ms  [A/m]
A=3.5e-12;      % exchange constant  [J/m]
u=pi*4e-7;      % permeability of vacuum  [H/m]
r=2.21e5;       % gyromagnetic ratio  [m/(s.A)]

j=1;

H=B*796;        % external field  [A/m] 
                % H=B/u0, u0=4*pi*e-7  [H/m]

damping=2e-4;   % damping
dH0=0.2*796;    % inhomogeneous linewidth
    
Wm=r*Ms*1e-9;   % [GHz]
Wh=r*H*1e-9;    % [GHz]
Le=sqrt(2*A/(u*Ms^2))*1e9;   % exchange length

weff=inf;       % effective width extract from Mumax3

k=1*pi/weff;
Ky=k;           %the effective wave number describing SW mode across the width direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%