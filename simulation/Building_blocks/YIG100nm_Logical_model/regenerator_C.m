function [out_signal] = regenerator_C(in_signal,model)

% This function describes the behavior of the regenerator (DC+amplifier)
% fot the carry bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad], delay]

out_signal = in_signal;

end

