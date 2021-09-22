% N-bit Ripple-carry adder (RCA_nbit)
full_adder

% Area
area_RCA_nbit = N * area_FA;  % [um^2]

    
% Propagation delay
tpd_RCA_nbit = (N-1) * tpd_FA_C + max(tpd_FA_C,tpd_FA_S);  % [ns]


% Energy
switch model
    case 'YIG 100nm'
        energy_RCA_nbit = 0.5*(N+1)*energy_sw + N*(energy_FA - 0.5*(3*energy_sw));  % [aJ]
    case 'YIG 30nm'
        energy_RCA_nbit = 0.5*(N+1)*energy_sw + N*(energy_FA - 0.5*(3*energy_sw));  % [aJ]
end

switch model
    case 'YIG 100nm'
        fprintf('(YIG - 100nm) %d-bit RCA: area = %d um^2, delay = %d ns, energy = %d aJ \n',N, area_RCA_nbit, tpd_RCA_nbit, energy_RCA_nbit)
    case 'YIG 30nm'
        fprintf('(YIG - 30nm) %d-bit RCA: area = %d um^2, delay = %d ns, energy = %d aJ \n',N, area_RCA_nbit, tpd_RCA_nbit, energy_RCA_nbit)
end
    