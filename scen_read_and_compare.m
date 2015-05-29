% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;
%input_name = 'data/impact';
% input_name1 = 'data/impact_bc';
% input_name2 = 'data/impact_org';

input_name1 = 'data/square3_org';
input_name2 = 'data/square3_bc';

ps1 = sph_scenario(input_name1);
ps2 = sph_scenario(input_name2);

%compare only every di-step
di = 1;

var_name = 'pj';
Neval   = 4000;
Nss = 2;
ps1.Neval = Neval;
ps2.Neval = Neval;
ps1.Nss = Nss;
ps2.Nss = Nss;


plotall = true;
ps1.plot_quantity = 'vpd';
ps2.plot_quantity = 'vpd';
ps1.fixaxes.p= 1e8 *[-2,2];
ps2.fixaxes.p= 1e8 *[-2,2];

%% create particle class
obj_particles1 = sph_particles(ps1);
obj_particles2 = sph_particles(ps2);

%
if plotall
    obj_particles1.IO.initialize();
    obj_particles2.IO.initialize();
end

errfig=figure();

data = h5info([input_name1,'.h5'],'/');

Omega_compare = [ps1.Omega(1,:);     %x
                 0,98e-3];%y 
             
% offset = 1e-2; %border particles will differ form the influence of out of border particles
% Omega_compare = ps2.geo.omega_geo +[offset,-offset;offset,-offset]; 

[x_eval,dx] = ps2.rectangle(Omega_compare,Neval);
Neval = size(x_eval,1);
Ncomp=Neval; %evaluation points to compare

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
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_name2,group);

    if plotall
       obj_particles1.IO.do(obj_particles1);
       obj_particles2.IO.do(obj_particles2);
    end

    dat_eval1 = obj_particles1.IO.eval (obj_particles1, var_name,x_eval);
    dat_eval2 = obj_particles2.IO.eval (obj_particles2, var_name,x_eval);
    
        %% evaluate differences: (with suppersampling)
%     dat_eval1 = obj_particles1.IO.eval_ss (obj_particles1, var_name);
%     dat_eval2 = obj_particles2.IO.eval_ss (obj_particles2, var_name);
%     
    figure(errfig);
    subplot(2,1,1)    
    %2d
%     plot(x_eval1,dat_eval1,...
%          x_eval2,dat_eval2);
    %3d
%     plot3(x_eval(:,1),x_eval(:,2),dat_eval1,'x',...
%           x_eval(:,1),x_eval(:,2),dat_eval2,'o');
    obj_particles1.IO.plot_trisurf(x_eval(:,1),x_eval(:,2),(dat_eval2-dat_eval1),2*max(dx));        
    colorbar 
    colormap jet
    view([1,-1,1])
    title([group,' ; ', var_name]);
    % time evolution of the l2-error
    subplot(2,1,2)
    t=str2double(group(2:end));
    abs_error = sum(((dat_eval1(1:Ncomp)-dat_eval2(1:Ncomp))).^2,1).^0.5; %l2-norm
    rel_error = abs_error / sum(((dat_eval1(1:Ncomp))).^2,1).^0.5;
    plot (t,abs_error/Ncomp,'x');
    hold on;
    title([group,'; ',num2str(Ncomp),' evaluation particles']);
    drawnow    
    
end
if plotall
    obj_particles1.IO.finalize();
    obj_particles2.IO.finalize();
end


