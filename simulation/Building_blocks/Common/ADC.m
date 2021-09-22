function [out_value] = ADC(in_signal,model_parameters)


%SW_parameters

N = size(in_signal);
N = N(1);

out_value = zeros(1,N);
for j=1:N
    if in_signal(j,1) < model_parameters.SW_amplitude/sqrt(3)
        out_value(j) = 0;
    elseif in_signal(j,1) > model_parameters.SW_amplitude/sqrt(3/2)
        out_value(j) = 1;
    else
        out_value(j) = 0.5;
        disp('ADC: forbidden state')
    end
end

end