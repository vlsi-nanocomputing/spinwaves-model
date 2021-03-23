function [out_value] = ADC(in_signal,model)


SW_parameters

N = size(in_signal);
N = N(1);

out_value = zeros(1,N);
for j=1:N
    if in_signal(j,1) < SW_amplitude/sqrt(3)
        out_value(j) = 0;
    else
        out_value(j) = 1;
    end
end

end