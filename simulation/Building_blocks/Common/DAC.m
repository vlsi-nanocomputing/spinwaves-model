function [out_signal] = DAC(in_value,model_parameters)
% Spin-wave generator (antenna)


%SW_parameters % script

N = size(in_value);
N = max(N);
out_signal = zeros(N,model_parameters.N_inf);
for j=1:N
    switch in_value(j)
        case 0
            out_signal(j,1) = 0;
        case 1
            out_signal(j,1) = model_parameters.SW_amplitude;
        otherwise
            disp('Error: check your digital inputs')
            out_signal(j,1) = inf;
    end
end
out_signal(:,2) = model_parameters.SW_frequency;
out_signal(:,3) = 0;
out_signal(:,4) = 0;

end

