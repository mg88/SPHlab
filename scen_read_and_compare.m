% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare two solutions

close all; clear; clc;

%input_name = 'data/impact';
input_name1 = 'data/impactISO_bc';
input_name2 = 'data/impactISO_org';

% input_name1 = 'data/square3_org';
% input_name2 = 'data/square3_bc';

%compare only every di-step
di = 1;
plot_solutions = false;
compare_solutions = false;
boundary_evolution = true;

plot_name  = 'p'; %what value shall I plot
plot_style = 'trisurf';
limaxes    = 1e8 *[-4,4];
movie_name = ''; %'movies/comp3.avi';  % empty means no video

var_name = 'pj'; %for evaluation
Neval   = 2000;
Nss = 2; %supersampling
Omega_compare = [-10e-3,30e-3;     %x
                 -60e-3,60e-3];%y 
             
x_evolution = [0.01, 0.03];             
    
%% generate objects:
ps1 = sph_scenario(input_name1);
ps2 = sph_scenario(input_name2);

%% write to object
ps1.Neval = Neval;
ps2.Neval = Neval;
ps1.Nss = Nss;
ps2.Nss = Nss;
ps1.plot_quantity = plot_name;
ps1.plot_style.p  = plot_style;
ps2.plot_style.p  = plot_style;
ps2.plot_quantity = plot_name;
ps1.fixaxes.p= limaxes;
ps2.fixaxes.p= limaxes;

%% create particle class
obj_particles1 = sph_particles(ps1);
obj_particles2 = sph_particles(ps2);
            
%% evaluation grid:
[x_eval,dx] = ps2.rectangle(Omega_compare,Neval);
Neval = size(x_eval,1); %real amount of evaluation points

%% what are the timesteps to read (the hdf5-groupname)
data = h5info([input_name1,'.h5'],'/');
i=1;
T=[];
while(i<size(data.Groups,1));    
    T(i) =  str2double(data.Groups(i).Name(2:end));
    i=i+1;
end
[~,I]=sort(T);
I = I(1:di:size(I,2));

%% create figure
errfig=figure();
nisubplot = plot_solutions + compare_solutions + boundary_evolution;
njsubplot = 2;

%% movie
if ~isempty(movie_name)
    set(gca,'DataAspectRatio',[1,1,1]);
    errfig.Units = 'normalized';
    errfig.Position = [0 0 1 1];
    vidObj  = VideoWriter(movie_name);
    open(vidObj);
end   

%% start reading input
nI = size(I,2);
T  = zeros(nI,1);
dat1_evolution = zeros(nI, 1);
dat2_evolution = zeros(nI, 1);

for i=1:nI;
    iplot = 1;
    %read data
    group = data.Groups(I(i)).Name;   %(time)
    t = str2double(group(2:end));
    T(i)=t;
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_name2,group);

    %% plot the solutions:
    if plot_solutions
        obj_particles1.IO.mfigure = errfig;
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        obj_particles1.IO.plot_data(obj_particles1,ps1.plot_quantity,false);
        obj_particles2.IO.mfigure = errfig;
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        obj_particles2.IO.plot_data(obj_particles2,ps2.plot_quantity,false);
    end

    %% evaluate two simulations on evaluation points and compare the results
    if compare_solutions
        dat_eval1 = obj_particles1.IO.eval (obj_particles1, var_name,x_eval);
        dat_eval2 = obj_particles2.IO.eval (obj_particles2, var_name,x_eval);


        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        %2d
    %     plot(x_eval1,dat_eval1,...
    %          x_eval2,dat_eval2);
        %3d
    %     plot3(x_eval(:,1),x_eval(:,2),dat_eval1,'x',...
    %           x_eval(:,1),x_eval(:,2),dat_eval2,'o');
    
        % plot difference
        obj_particles1.IO.plot_trisurf(x_eval(:,1),x_eval(:,2),(dat_eval2-dat_eval1),2*max(dx));        
        colorbar 
        caxis(limaxes)
        colormap jet
        view([1,-1,1])
        title([group,' ; ', var_name]);
        
        % time evolution of the l2-error
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        abs_error = sum(((dat_eval1(1:Neval)-dat_eval2(1:Neval))).^2,1).^0.5; %l2-norm
        rel_error = abs_error / sum(((dat_eval1(1:Neval))).^2,1).^0.5;
        plot (t,abs_error/Neval,'x');
        hold on;
        title([group,'; ',num2str(Neval),' evaluation particles']);
    end
    
    if boundary_evolution
        dat1_evolution(i,:) = obj_particles1.IO.eval (obj_particles1, var_name,x_evolution);
        dat2_evolution(i,:) = obj_particles2.IO.eval (obj_particles2, var_name,x_evolution);      
    end
    
    drawnow        
    %% movie:
    if ~isempty(movie_name)
        currFrame = getframe(errfig);
        writeVideo(vidObj,currFrame);
    end            
end
if boundary_evolution
    figure 
    plot(T,dat1_evolution,'b-',...
         T,dat2_evolution,'r-');   
    title('evolution at point x-evolution');
    legend(input_name1,input_name2);
    xlabel('t');
    ylabel(var_name);
    keyboard
end

if ~isempty(movie_name)
  close(vidObj);
  disp (['movie saved as ',movie_name]);
end   
