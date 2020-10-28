function [output_power] = amplifier(input_power)




% In order to compensate the energy splitting in the DC1 and parasitic
% reflections and spin-wave damping in the waveguides, the correct gain of
% the amplifier is higher than 2. In these simplified models, the
% waveguide is lossless, so the gain=2 that compensate only the energy
% splitting
output_power = input_power * 2;

end

