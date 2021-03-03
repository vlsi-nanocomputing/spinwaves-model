function [out_signal] = amplifier(in_signal,gain)

% This function describes the relation between the input and output of a
% VCMA amplifier.
% Using the 'optimal_gain_ver2.m' script, you can see what is the effect of
% the gain on the '0' and '1' logic. You can find a desired gain.

% The input and output variables are vectors, and they are composed in
% the following way:
% [amplitude(dimensionless), frequency [GHz], phase [rad]]

out_signal = in_signal;
out_signal(1) = in_signal(1)*sqrt(gain);
end

