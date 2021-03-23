function [norm_powers] = normalization(in_signals,model)

% The function receives the signals amplitude, and it gives the normalized power.

SW_parameters % script

norm_powers = (in_signals/SW_amplitude).^2*100;  %  (%)

