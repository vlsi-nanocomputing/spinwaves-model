function [out_signal] = regenerator_S(in_signal,model)

% This function describes the behavior of the regenerator (DC+amplifier)
% for the sum bit output.
% This block regenerates the correct SW amplitude according to the logic
% value.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad]]


out_signal = amplifier(in_signal,4);

end

