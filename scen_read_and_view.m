% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;
%input_name = 'data/impact';
input_name1 = 'data/riemann2d_inf';

ps1 = sph_scenario(input_name1);

ps1.plot_quantity = 'vpd';


%% create particle class
obj_particles1 = sph_particles(ps1);

%
obj_particles1.IO.initialize();


data = h5info([input_name1,'.h5'],'/');

i=1;
while(i<size(data.Groups,1));    
    %read data
    group = data.Groups(i).Name;   %(time)
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);

    obj_particles1.IO.do(obj_particles1);

    %% next step
    i=i+1;
end
obj_particles1.IO.finalize();


