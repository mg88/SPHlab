% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;

input_name = 'impactMG_monolithic';
in_namedir = ['data/',input_name];
% input_name1 = 'data/square3_bc';

% input_name1 = '../../../../Desktop/out-hvi';


ps1 = sph_scenario(in_namedir);
ps1.plotconfig.figurename = input_name;
ps1.plotconfig.latexplot = true;
% ps1.plot_quantity = 'xp';
% ps1.plot_style.x = 'ringcloud';
ps1.plot_style.e = 'patches';
%


% intervall of plotting (every di step)
di = 2;

% ps1.eta = 1.2;
% ps1.Omega = [ -0.1, 0.4;
%             -0.3, 0.3];
% ps1.plot_dt = 1e-6;  


%% create particle class
obj_particles1 = sph_particles(ps1);

%
obj_particles1.IO.initialize();


data = h5info([in_namedir,'.h5'],'/');

i=1;
T=[];
while(i<size(data.Groups,1));    
    T(i) =  str2double(data.Groups(i).Name(2:end));
    i=i+1;
end
[~,I]=sort(T);

I = I(1:di:size(I,2));

for i=I;    
    %read data
    group = data.Groups(i).Name;   %(time)
    obj_particles1.IO.read_hdf5(obj_particles1,in_namedir,group);    
    obj_particles1.IO.do(obj_particles1,true,false);
%     obj_part#icles1.IO.plot_data(obj_particles1,'p'); 
    obj_particles1.IO.savefigure(group(2:end),'','png')
end
obj_particles1.IO.finalize();


