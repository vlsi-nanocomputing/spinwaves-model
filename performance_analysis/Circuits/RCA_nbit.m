% N-bit Ripple-carry adder (RCA_nbit)
full_adder

% N = 32; % parallelism

% Area
area_RCA_nbit = N * area_FA;  % [um^2]

    
% Propagation delay
tpd_RCA_nbit = (N-1) * tpd_FA_C + max(tpd_FA_C,tpd_FA_S);  % [ns]


% Energy
if model == 1  % CMOS
    Ps_RCA_nbit = N*ps_FA; % [W]
    energy_RCA_nbit = Ps_RCA_nbit * tpd_RCA_nbit * 1e9 + N * pdyn_FA * 1e18; % [aJ]
else
    energy_RCA_nbit = 0.5*(N+1)*energy_sw + N*(energy_FA - 0.5*(3*energy_sw));  % [aJ]
end    


if model == 1
    fprintf('(CMOS - 45nm) %d-bit RCA: area = %d um^2, delay = %d ns, energy = %d aJ \n',N, area_RCA_nbit, tpd_RCA_nbit, energy_RCA_nbit)
elseif model == 2
    fprintf('(YIG - 100nm) %d-bit RCA: area = %d um^2, delay = %d ns, energy = %d aJ \n',N, area_RCA_nbit, tpd_RCA_nbit, energy_RCA_nbit)
else
    fprintf('(YIG - 30nm) %d-bit RCA: area = %d um^2, delay = %d ns, energy = %d aJ \n',N, area_RCA_nbit, tpd_RCA_nbit, energy_RCA_nbit)
end
    