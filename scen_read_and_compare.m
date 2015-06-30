% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare two solutions

% close all; clear; clc;

%input_name = 'data/impact';
input_name1 = 'data/impactISO_org2';
input_name2 = 'data/impactISO_bc2';

% input_name1 = 'data/square3_org';
% input_name2 = 'data/square3_bc';

%compare only every di-step
di = 1;
Tstart = 0.e-4;%0.68e-4;
Tend = 0.55e-4;
compare_position   = 0;
plot_solutions     = 0;
plot_difference    = 0;
l2error_evolution  = 0;
boundary_evolution = 1;

% zoom in:
zoom_in = true;
x_zoom = [0.001,0.03];
r_zoom = 0.005;

% compare position
extra_magnify = false;
x_magnify = [0.01,0.03];%abs(x(Aistar(1)));
r_magnify = 0.005;
            
% compare solution
plot_name  = 'p'; %what value shall I plot
plot_style = 'trisurf';
limaxes    = 1e8 *[-8,8];
movie_name = '';%'test'; %'movies/comp3.avi';  % empty means no video

var_name = 'pj'; %for evaluation
Neval   = 2000;
Nss = 1; %supersampling
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
ps1.plot_style.p  = plot_style;
ps2.plot_style.p  = plot_style;
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

% take only the time smaller Tend and greater than Tstart
I=I(logical((T(I)>Tstart) .* (T(I)<Tend)));

%% create figures
% figure for position comparison  (in one plot)
position_fig = figure();
if extra_magnify
    %get data from original axes
    ax1=gca;
    set(ax1,'units','pixels');
    pos=get(ax1,'Position');

    %plot second axis
    dx=pos(3)/12;
    dy=pos(4)/12;
    b = pos(3)/2-dx;
    h = pos(4)/2-dy;
    x00=pos(1)+dx;
    y00=pos(2)+pos(4)/2 + dy;

    ax2=axes();
    set(ax2,'units','pixels');
    set(ax2,'position',[x00,y00,min(b,h),min(b,h)]);
    axis equal
    box on
    %zoom in        
    xlim ([x_magnify(1)-r_magnify,x_magnify(1)+r_magnify]);
    ylim ([x_magnify(2)-r_magnify,x_magnify(2)+r_magnify]);
end

% figure for variable comparison (in subplots) + errorplot
nisubplot = 1;
njsubplot = 2*plot_solutions + plot_difference;
if njsubplot >0
    errfig = figure();
end

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
l2err_evolution = zeros(nI, 1);


for i=1:nI;
    iplot = 1;
    %read data
    group = data.Groups(I(i)).Name;   %(time)
    t = str2double(group(2:end));
    T(i)=t;
    obj_particles1.IO.read_hdf5(obj_particles1,input_name1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_name2,group);

    %compare the positions of the particles
    if compare_position
        if extra_magnify
            axes(ax1);
        end
        figure(position_fig)
        obj_particles1.IO.mfigure = position_fig;
        obj_particles2.IO.mfigure = position_fig;
        %plot position (with clearing axes) -> o-marker
        obj_particles1.IO.plot_data(obj_particles1,'x',false,true);
        %plot position (without clearing axes) -> x-marker
        obj_particles2.IO.plot_data(obj_particles2,'x',false,false);
        set(gca,'DataAspectRatio',[1,1,1]);
        
        if zoom_in
            xlim ([x_zoom(1)-r_zoom,x_zoom(1)+r_zoom]);
            ylim ([x_zoom(2)-r_zoom,x_zoom(2)+r_zoom]);
        end
        
        if extra_magnify
            %second axis (magnification)
            axes(ax2)
            cla
            copyobj(allchild(ax1),ax2);
            %draw rectangle (not earlier otherwise it would be copied too.
            axes(ax1);
            rectangle('Position',[x_magnify(1)-r_magnify,x_magnify(2)-r_magnify,x_magnify(1)+r_magnify,x_magnify(1)+r_magnify])
            axes(ax2);
        end
    end
    
    
    %% plot the solutions:
    if plot_solutions
        figure(errfig)
        obj_particles1.IO.mfigure = errfig;
        obj_particles2.IO.mfigure = errfig;
        %plot 1
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        obj_particles1.IO.plot_data(obj_particles1,plot_name,false);                
        if zoom_in
            xlim ([x_zoom(1)-r_zoom,x_zoom(1)+r_zoom]);
            ylim ([x_zoom(2)-r_zoom,x_zoom(2)+r_zoom]);
        end
        %plot 2 (original)
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        obj_particles2.IO.plot_data(obj_particles2,plot_name,false);
        
         % ----------
        %draw same boundary line as above
        hold on;
        obj_particles1.IO.draw_geometry(obj_particles1)

%         %boundary lines:
%         hold on;
%         a = axis; %now with Omega in consideration
%         z=0;
%         for boun = obj_particles1.bc
%            e = [0 -1; 1 0]*boun.outer_normal'; %rotate 90degrees
%            p1= boun.bp +1000*e'; %from "minus infinity to infinity"
%            p2= boun.bp -1000*e';
%            plot3([p1(1),p2(1)],[p1(2),p2(2)],[z,z],'kx:');
%         end
%         axis(a);
        % ----------
        
        if zoom_in
            xlim ([x_zoom(1)-r_zoom,x_zoom(1)+r_zoom]);
            ylim ([x_zoom(2)-r_zoom,x_zoom(2)+r_zoom]);
        end
    end

    %% evaluate two simulations on evaluation points and compare the results
    if plot_difference || l2error_evolution
        dat_eval1 = obj_particles1.IO.eval (obj_particles1, var_name,x_eval);
        dat_eval2 = obj_particles2.IO.eval (obj_particles2, var_name,x_eval);
    end
    
    if plot_difference
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        cla
        %2d
    %     plot(x_eval1,dat_eval1,...
    %          x_eval2,dat_eval2);
        %3d
    %     plot3(x_eval(:,1),x_eval(:,2),dat_eval1,'x',...
    %           x_eval(:,1),x_eval(:,2),dat_eval2,'o');
    
        % plot difference           
        obj_particles1.IO.plot_trisurf(x_eval(:,1),x_eval(:,2),(dat_eval2-dat_eval1),2*max(dx));        
        hold on;
        %add Omega and boundary lines
        obj_particles1.IO.draw_geometry(obj_particles1)
        set(gca,'gridlinestyle','none');
        caxis(limaxes)
        view (2)
        daspect([1,1,1])
        colorbar 
        colormap jet
        title(['error at ',group,' ; ', var_name]);
    end
    
    %save l2error
    if l2error_evolution
        % time evolution of the l2-error
        abs_error = sum(((dat_eval1(1:Neval)-dat_eval2(1:Neval))).^2,1).^0.5; %l2-norm
%         rel_error = abs_error / sum(((dat_eval1(1:Neval))).^2,1).^0.5;
        l2err_evolution(i) = abs_error/Neval;
    end
    
    if boundary_evolution
        dat1_evolution(i,:) = obj_particles1.IO.eval (obj_particles1, var_name,x_evolution);
        dat2_evolution(i,:) = obj_particles2.IO.eval (obj_particles2, var_name,x_evolution);      
    end
    
    drawnow        
    %% movie:
    if ~isempty(movie_name)
         currFrame = getframe(errfig);
%         writeVideo(vidObj,currFrame);
          im = frame2im(currFrame);
          [imind,cm] = rgb2ind(im,256);
          if i == 1;
              imwrite(imind,cm,movie_name,'gif', 'Loopcount',inf);
          else
              imwrite(imind,cm,movie_name,'gif','WriteMode','append');
          end
        
    end            
end

%%
if l2error_evolution
    figure
    plot (T,l2err_evolution,'xb-');
    title([group,'; ',num2str(Neval),' evaluation particles']);   
end

%%
if boundary_evolution
    figure 
    plot(T,dat1_evolution,'b-',...
         T,dat2_evolution,'r-');   
    title('evolution at point x');
    legend(input_name1,input_name2);
    xlabel('t');
    ylabel(var_name);
end

if ~isempty(movie_name)
  close(vidObj);
  disp (['movie saved as ',movie_name]);
end   
