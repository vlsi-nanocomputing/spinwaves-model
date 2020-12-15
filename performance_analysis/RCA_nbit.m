% N-bit Ripple-carry adder (RCA_nbit)
full_adder

% N = 2; % parallelism

% Area
A_RCA_nbit = N * A_FA;  % [um^2]

% Energy
E_RCA_nbit = (N+0.5) * E_a_b + N*3*E_amp_regS;  % [aJ]

% Propagation delay
tpd_RCA_nbit = (N-1) * tpd_FA_C + max(tpd_FA_C,tpd_FA_S);  % [ns]



