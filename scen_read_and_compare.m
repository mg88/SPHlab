% SPH for HVI  - Markus Ganser - TU/e - 2015
% read and compare two solutions

% close all; clear; clc;

% input names:
input_name1 = 'impact_laminate_test';
input_name2 = 'impact_laminate_test';
%--------------------


input_namedir1 = ['data/',input_name1];
input_namedir2 = ['data/',input_name2];


%compare only every di-step
di = 1;
Tstart = 0;%1.9e-5;%0.68e-4;
Tend = 5.e-6;
compare_position   = 0;
plot_solutions     = 1;
plot_difference    = 0;
l2error_evolution  = 0;
boundary_evolution = 0;

% zoom in:
zoom_in = true;
% x_zoom = [-0.5e-3,12e-3]; %co
% r_zoom = 1.7e-3;
x_zoom = [3e-3,0];
r_zoom = 22e-3;

% compare position
extra_magnify = false;
x_magnify = [-0.5e-3,11e-3];%abs(x(Aistar(1)));
r_magnify = 1.5e-3;
            
% compare solution (show side by side)
plot_name  = 'x'; %what value shall I plot
plot_style = 'patches';
limaxes    = '';%1e7 *[-8,8];
movie_name = '';%'movies/impact_comp2.mp4';%'test'; %'movies/comp3.avi';  % empty means no video

var_name = 'pj'; %for evaluation
Neval   = 2000;
Nss = 1; %supersampling
Omega_compare = [-5e-3,20e-3;     %x
                 -20e-3,20e-3];%y 
             
% boundary evolution: evaluation point:
x_evolution = [3e-3,8e-3];%[2.5e-3, 15e-3];             
    
%% generate objects:
ps1 = SPHlab_scenario(input_namedir1);
ps2 = SPHlab_scenario(input_namedir2);
%
ps1.Omega = Omega_compare;
ps2.Omega = Omega_compare;
ps1.plotconfig.figuresize = []; %empty: make default
ps1.plotconfig.transpose = []; %empty: make default
ps1.plotconfig.latexplot = true; %empty: make default
ps2.plotconfig.figuresize = []; %empty: make default
ps2.plotconfig.transpose = []; %empty: make default
ps2.plotconfig.latexplot = true; %empty: make default

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
obj_particles1 = SPHlab_particles(ps1);
obj_particles2 = SPHlab_particles(ps2);
            
%% evaluation grid:
[x_eval,~] = ps2.rectangle(Omega_compare,Neval);
Neval = size(x_eval,1); %real amount of evaluation points

%% what are the timesteps to read (the hdf5-groupname)
data = h5info([input_namedir1,'.h5'],'/');
data2 = h5info([input_namedir2,'.h5'],'/');
i=1;
T=[];
while(i<size(data.Groups,1));    
    T(i) =  str2double(data.Groups(i).Name(2:end));
    i=i+1;
end
i=1;
T2=[];
while(i<size(data2.Groups,1));    
    T2(i) =  str2double(data2.Groups(i).Name(2:end));
    i=i+1;
end
[~,I]=sort(T);
I = I(1:di:size(I,2));
[~,I2]=sort(T2);
I2 = I2(1:di:size(I2,2));
if size(T,2)~=size(T2,2)
    warning('take care, date might be not defined on the same timestep')
end

% take only the time smaller Tend and greater than Tstart
I=I(logical((T(I)>Tstart) .* (T(I)<Tend)));
I2=I2(logical((T2(I2)>Tstart) .* (T2(I2)<Tend)));

%% create figures
% figure for position comparison  (in one plot)
position_fig = figure('Color', [1 1 1]);
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
    errfig = figure;%('Color',[1 1 1]);
end

%% movie
if ~isempty(movie_name)
    set(gca,'DataAspectRatio',[1,1,1]);
    errfig.Units = 'normalized';
    errfig.Position = [0.1 0.1 0.9 0.9];
    vidObj  = VideoWriter(movie_name,'MPEG-4');
    vidObj.FrameRate = 5;
    open(vidObj);
end   

%% start reading input
nI = min(size(I,2),size(I2,2));
T  = zeros(nI,1);
dat1_evolution = zeros(nI, 1);
dat2_evolution = zeros(nI, 1);
l2err_evolution = zeros(nI, 1);


for i=1:nI;
    iplot = 1;
    %read data
    group = data.Groups(I(i)).Name;   %(time)
    group2 = data2.Groups(I2(i)).Name;   %(time)

    t = str2double(group(2:end));
    T(i)=t;
    obj_particles1.IO.read_hdf5(obj_particles1,input_namedir1,group);
    obj_particles2.IO.read_hdf5(obj_particles2,input_namedir2,group2);

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
          % ----------
        %draw same boundary line as above
        hold on;
        obj_particles1.IO.draw_geometry(obj_particles2)
        
        if zoom_in
            xlim ([x_zoom(1)-r_zoom,x_zoom(1)+r_zoom]);
            ylim ([x_zoom(2)-r_zoom,x_zoom(2)+r_zoom]);
        end
        %plot 2 (original)
        subplot(nisubplot,njsubplot,iplot); iplot=iplot+1;
        obj_particles2.IO.plot_data(obj_particles2,plot_name,false);
        


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
%         % ----------
%         
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
        obj_particles1.IO.draw_geometry(obj_particles2)
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
         writeVideo(vidObj,currFrame);
%           im = frame2im(currFrame);
%           [imind,cm] = rgb2ind(im,256);
%           if i == 1;
%               imwrite(imind,cm,movie_name,'gif', 'Loopcount',inf);
%           else
%               imwrite(imind,cm,movie_name,'gif','WriteMode','append');
%           end
        
    end            
end

% sompe layout properties: (plot in latex style)
%---------------- 
FontSize=15;
Sfont= struct('FontUnits','points',...
'FontSize',FontSize,...
'FontName','Times');
Sint = struct('interpreter','latex');

%%
if l2error_evolution
     figure('Color',[1 1 1]); 
    plot (T,l2err_evolution,'xb-');
    title([group,'; ',num2str(Neval),' evaluation particles']);   
end

%%
if boundary_evolution
    figure('Color',[1 1 1]); 
    plot(T,dat1_evolution,'b-',...
         T,dat2_evolution,'r-');   
%     title('evolution at point x');
    h=legend(input_namedir1,input_namedir2);
    set(h,Sfont,Sint);

    h1=xlabel('$t$');
    h2=ylabel(var_name);
    set(h1,Sfont,Sint);
    set(h2,Sfont,Sint);
    set(gca,Sfont);

end

if ~isempty(movie_name)
  close(vidObj);
  disp (['movie saved as ',movie_name]);
end   
