function [out_signal] = amplifier(in_signal,gain)

% This function describes the relation between the input and output of a
% VCMA amplifier.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad], delay [ns]]

N = size(in_signal);
N = N(2);

if N == 1 % N=1 is the logical model
	out_signal = in_signal * gain;
else
	out_signal = in_signal;
	out_signal(1) = in_signal(1)*sqrt(gain);

end

end

