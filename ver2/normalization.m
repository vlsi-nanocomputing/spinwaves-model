function [norm_powers] = normalization(in_signals)

% The function receives the signals amplitude, and it gives the normalized power.
SW_amplitude = 0.153; 

norm_powers = (in_signals/SW_amplitude).^2*100;

