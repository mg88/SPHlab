% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;


%% input name:
input_name = 'impact_laminate_test';


in_namedir = ['data/',input_name];


ps1 = SPHlab_scenario(in_namedir);
ps1.plotconfig.figurename = input_name;
ps1.plotconfig.latexplot = true;
ps1.plot_quantity = 'x';
ps1.fixaxes.p = 5*[-1e8,1e8];
ps1.fixaxes.e = 1*[-1e5,1e5];
ps1.fixaxes.d = [1000,2900];

ps1.plot_style.d = 'patches';

% ps1.plot_style.x = 'ringcloud';
ps1.plotconfig.figuresize = []; %empty: make default
ps1.plotconfig.transpose = []; %empty: make default


% ps1.plot_style.x = 'ringcloud';
% ps1.plot_style.e = 'patches';
ps1.fixaxes.p = [-5e8,5e8];

%movie:
ps1.save_as_movie = false;
ps1.movie_name = input_name;

%figure:
ps1.save_as_figure = false;
ps1.figure_format = 'eps';

% intervall of plotting (every di step)
di = 2;

% ps1.eta = 1.2;
% ps1.Omega = [ -0.05, 0.1;
%             -0.25, 0.25];
% ps1.plot_dt = 1e-6;  


%% create particle class
obj_particles1 = SPHlab_particles(ps1);

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
%     obj_particles1.IO.savefigure(group(2:end),'','png')
    disp(group)   
%         xlim([-0.015,0.02]);
    
    
end
obj_particles1.IO.finalize();


