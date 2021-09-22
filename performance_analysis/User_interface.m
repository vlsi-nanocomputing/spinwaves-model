% User interface

model = 'YIG 30nm';

switch model 
    case 'YIG 100nm'
        model_path = 'YIG/YIG_100nm';
        addpath('YIG')
    case 'YIG 30nm'
        model_path = 'YIG/YIG_30nm';
        addpath('YIG')
end
addpath(model_path)

N = input('Insert the parallelism of the RCA \n');

addpath('Circuits')
RCA_nbit




