% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare solutions

close all; clear; clc;
%input_name = 'data/impact';
input_name1 = 'data/impactISO_bc';
input_name2 = 'data/impactISO_org';

% input_name1 = 'data/square3_org';
% input_name2 = 'data/square3_bc';

ps1 = sph_scenario(input_name1);
ps2 = sph_scenario(input_name2);

%compare only every di-step
di = 2;

var_name = 'pj';
plot_name = 'p';
movie_name = '';%movies/comp2.avi';
Neval   = 2000;
Nss = 2;
ps1.Neval = Neval;
ps2.Neval = Neval;
ps1.Nss = Nss;
ps2.Nss = Nss;

plotall = false;
ps1.plot_quantity = plot_name;
ps1.plot_style.p='patches';
ps2.plot_style.p='patches';
ps2.plot_quantity = plot_name;
limaxes=1e8 *[-4,4];
ps1.fixaxes.p= limaxes;
ps2.fixaxes.p= limaxes;

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

Omega_compare = [-10e-3,30e-3;     %x
                 -60e-3,60e-3];%y 
             
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

%movie
if ~isempty(movie_name)
    set(gca,'DataAspectRatio',[1,1,1]);
    errfig.Units = 'normalized';
    errfig.Position = [0 0 1 1];
    vidObj  = VideoWriter(movie_name);
    open(vidObj);
end   

for i=I;
    %read data
    group = data.Groups(i).Name;   %(time)
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_name2,group);

    %plot solutions
    obj_particles1.IO.mfigure = errfig;
    subplot(2,2,1);
    obj_particles1.IO.plot_data(obj_particles1,ps1.plot_quantity);
    obj_particles2.IO.mfigure = errfig;
    subplot(2,2,2);
    obj_particles2.IO.plot_data(obj_particles2,ps2.plot_quantity);

    dat_eval1 = obj_particles1.IO.eval (obj_particles1, var_name,x_eval);
    dat_eval2 = obj_particles2.IO.eval (obj_particles2, var_name,x_eval);
    
 
%     figure(errfig);
    subplot(2,2,3)    
    %2d
%     plot(x_eval1,dat_eval1,...
%          x_eval2,dat_eval2);
    %3d
%     plot3(x_eval(:,1),x_eval(:,2),dat_eval1,'x',...
%           x_eval(:,1),x_eval(:,2),dat_eval2,'o');
    obj_particles1.IO.plot_trisurf(x_eval(:,1),x_eval(:,2),(dat_eval2-dat_eval1),2*max(dx));        
    colorbar 
    caxis(limaxes)
    colormap jet
    view([1,-1,1])
    title([group,' ; ', var_name]);
    % time evolution of the l2-error
    subplot(2,2,4)
    t=str2double(group(2:end));
    abs_error = sum(((dat_eval1(1:Ncomp)-dat_eval2(1:Ncomp))).^2,1).^0.5; %l2-norm
    rel_error = abs_error / sum(((dat_eval1(1:Ncomp))).^2,1).^0.5;
    plot (t,abs_error/Ncomp,'x');
    hold on;
    title([group,'; ',num2str(Ncomp),' evaluation particles']);
    drawnow    
    
    %movie:
    if ~isempty(movie_name)
        currFrame = getframe(errfig);
        writeVideo(vidObj,currFrame);
    end    
    
    
end
if plotall
    obj_particles1.IO.finalize();
    obj_particles2.IO.finalize();
end

if ~isempty(movie_name)
  close(vidObj);
  disp (['movie saved as ',movie_name]);
end   
