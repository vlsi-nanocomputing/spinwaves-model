function [out_signal] = DAC(in_value)

N = size(in_value);
N = max(N);
out_signal = zeros(N,3);
for j=1:N
    switch in_value(j)
        case 0
            out_signal(j,1) = 0;
        case 1
            out_signal(j,1) = 0.153;
        otherwise
            display("Error: check your digital inputs")
            out_signal(j,1) = inf;
    end
end
out_signal(:,2) = 2.282;
out_signal(:,3) = 0;


end

