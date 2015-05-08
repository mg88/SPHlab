% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;
%input_name = 'data/impact';
input_name1 = 'data/riemann2d_inf';
input_name2 = 'data/riemann2d_boun'; %_mass

ps1 = sph_scenario(input_name1);
ps2 = sph_scenario(input_name2);

var_name = 'pj';
Neval   = 500;
Nss = 2;
ps1.Neval = Neval;
ps2.Neval = Neval;
ps1.Nss = Nss;
ps2.Nss = Nss;

Ncomp=Neval/2; %evaluation points to compare

plotall = false;
ps1.plot_quantity = 'vpd';
ps2.plot_quantity = 'vpd';


%% create particle class
obj_particles1 = sph_particles(ps1);
obj_particles2 = sph_particles(ps2);

%
if plotall
    obj_particles1.IO.initialize();
    obj_particles2.IO.initialize();
end


data = h5info([input_name1,'.h5'],'/');

x_eval1  = obj_particles1.IO.x_eval;
x_eval2  = obj_particles2.IO.x_eval;

i=1;
while(i<size(data.Groups,1));    
    %read data
    group = data.Groups(i).Name;   %(time)
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_name2,group);

    if plotall
       obj_particles1.IO.do(obj_particles1);
       obj_particles2.IO.do(obj_particles2);
    end

    dat_eval1 = obj_particles1.IO.eval (obj_particles1, var_name,x_eval1);
    dat_eval2 = obj_particles2.IO.eval (obj_particles2, var_name,x_eval2);
    
        %% evaluate differences: (with suppersampling)
%     dat_eval1 = obj_particles1.IO.eval_ss (obj_particles1, var_name);
%     dat_eval2 = obj_particles2.IO.eval_ss (obj_particles2, var_name);
%     
    
    subplot(2,1,1)    
    %2d
%     plot(x_eval1,dat_eval1,...
%          x_eval2,dat_eval2);
    %3d
    plot3(x_eval1(:,1),x_eval1(:,2),dat_eval1,'x',...
          x_eval2(:,1),x_eval2(:,2),dat_eval2,'o');
   view([0,-1,0])
    title([group,' ; ', var_name]);
    % time evolution of the l2-error
    subplot(2,1,2)
    t=str2double(group(2:end));
    abs_error = sum(((dat_eval1(1:Ncomp)-dat_eval2(1:Ncomp))).^2,1).^0.5; %l2-norm
    rel_error = abs_error / sum(((dat_eval1(1:Ncomp))).^2,1).^0.5;
    plot (t,abs_error/Ncomp,'x');
    hold on;
    title([group,'; first ',num2str(Ncomp),' particles']);
    drawnow    
    
    %% next step
    i=i+1;
end
if plotall
    obj_particles1.IO.finalize();
    obj_particles2.IO.finalize();
end


