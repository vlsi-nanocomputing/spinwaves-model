function [out_signal] =phase_shifter(in_signal,phase)

% This function describes the relation between the input and output of a
% phase shifter.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad]]

out_signal = in_signal;
out_signal(3) = in_signal(3) + phase;

end

