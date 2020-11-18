function [out_signal] = DAC(in_value)

switch in_value
    case 0
        out_signal(1) = 0;
    case 1
        out_signal(1) = 0.153;
    otherwise
        display("Error: check your digital inputs")
end
out_signal(2) = 2.282;
out_signal(3) = 0;

end

