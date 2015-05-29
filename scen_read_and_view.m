% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;

input_name1 = 'data/impact_bc';
% input_name1 = 'data/square3_bc';

% input_name1 = '../../../../Desktop/out-hvi';


ps1 = sph_scenario(input_name1);

ps1.plot_quantity = 'xvpde';
ps1.plot_style.p = 'trisurf';
ps1.plot_style.e = 'trisurf';
%
% ps1.eta = 1.2;
% ps1.Omega = [ -0.1, 0.4;
%             -0.3, 0.3];
% ps1.plot_dt = 1e-6;  


%% create particle class
obj_particles1 = sph_particles(ps1);

%
obj_particles1.IO.initialize();


data = h5info([input_name1,'.h5'],'/');

i=1;
T=[];
while(i<size(data.Groups,1));    
    T(i) =  str2double(data.Groups(i).Name(2:end));
    i=i+1;
end
[~,I]=sort(T);

di = 1;
I = I(1:di:size(I,2));

for i=I;    
    %read data
    group = data.Groups(i).Name;   %(time)
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);    
    obj_particles1.IO.do(obj_particles1);
%     obj_part#icles1.IO.plot_data(obj_particles1,'p'); 

end
obj_particles1.IO.finalize();


