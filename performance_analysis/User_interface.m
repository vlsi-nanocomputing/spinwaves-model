% User interface

model=2;
while model~=1 && model~=2 && model~=3 && model~=4 
    model = input('Choose one technology from the following list:\n 1) CMOS (45nm) \n 2) YIG (100nm) \n 3) YIG (30nm) \n 4) QUIT \n');
end
if model~=4

switch model 
    case 1
        model_path = 'CMOS_45nm';
    case 2
        model_path = 'YIG/YIG_100nm';
        addpath('YIG')
    case 3
        model_path = 'YIG/YIG_30nm';
        addpath('YIG')
end
addpath(model_path)

N = input('Insert the parallelism of the RCA \n');

addpath('Circuits')
RCA_nbit





end


