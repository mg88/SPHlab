classdef sph_IO < handle
    %IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %input
        
        %output
        dim
        
        %plot
        mfigure
        plotstyle
        plot_dt
        
        %movie
        save_as_movie
        movie_name
        vidObj
        
        %
        obj_particles
    end
    
    methods
        
%         %% constructor
        function obj = sph_IO(obj_scen)
              obj.save_as_movie = obj_scen.save_as_movie;
              obj.movie_name    = get_movie_name(obj_scen);
              obj.plot_dt       = obj_scen.plot_dt;
              obj.plotstyle     = obj_scen.plotstyle;
        end
        %%
        function initialize(obj)
            %% output
            disp(obj.dim)
            %% plot
            
            
            %% movie
            obj.mfigure = figure;
            if obj.save_as_movie
                set(gca,'DataAspectRatio',[1,1,1]);
                obj.vidObj     = VideoWriter(obj.movie_name);
                open(obj.vidObj);
            end            
            
        end
        %%
        function do (obj,t, obj_particles)
            %% write output
            obj.obj_particles= obj_particles;
            %%
            %% plotting
            if (mod(t,obj.plot_dt) < obj_particles.dt) %ToDo: good for variable timestepping?
                plot_data (obj,t);
                if obj.save_as_movie
                    currFrame = getframe(obj.mfigure);
                    writeVideo(obj.vidObj,currFrame);
                end
            end   
        end       
        %%
        function finalize (obj)
            %% close output
            
            %% movie
            if obj.save_as_movie
              close(obj.vidObj);
              disp (['movie saved as ',obj.movie_name]);
            end           
        end
        
        %% %%%  Plotting functions %%% %%
        function plot_data(obj,t)
           figure(obj.mfigure);
           if obj.obj_particles.dim == 1               
               plot_data1D(obj,t);
           elseif obj.obj_particles.dim == 2
               if strcmp(obj.plotstyle,'scatter')
                    plot_data2D(obj,t);
               elseif strcmp (obj.plotstyle,'trisurf')
                    plot_trisurf(obj,t)
               elseif strcmp (obj.plotstyle,'patches')
                    plot_patches(obj,t)
               else
                   error([obj.plotstyle, '- plotstyle is not supported']);
               end
           end 
        end          
        %%
        function plot_data1D(obj,time)
            data = obj.obj_particles;
            nplot = length(obj.plotstyle);
            iplot = 1;
            if ~isempty(strfind(obj.plotstyle,'x'));
                subplot(nplot,1,iplot)
                plot(data.Xj(data.Iin),0,'bo',data.Xj(data.Iboun),zeros(size(data.Iboun,1),1),'ko');
                title(['position, t=',num2str(time),' N= ',num2str(data.N)])
                xlim([data.Omega(1) data.Omega(2)]);
                iplot=iplot+1;     
            end
            if ~isempty(strfind(obj.plotstyle,'v'));
                iplot = create_a_supplot_1d(obj,data.vj,'velocity',time,nplot,iplot);
            end
            if ~isempty(strfind(obj.plotstyle,'p'));
                iplot = create_a_supplot_1d(obj,data.pj,'pressure',time,nplot,iplot);
            end
            if ~isempty(strfind(obj.plotstyle,'d'));
                iplot = create_a_supplot_1d(obj,data.rhoj,'density',time,nplot,iplot);
            end
            if ~isempty(strfind(obj.plotstyle,'f'));     
                subplot(nplot,1,iplot)
                bar(data.Xj,[data.F_int,data.F_diss,data.F_diss_art,data.F_ST]); title('forces')
                xlim([data.Omega(1) data.Omega(2)]);
                legend('int','diss','diss_art','st');
            end
            drawnow
        end
        %%
        function iplot = create_a_supplot_1d(obj,y,title_name,time,nplot,iplot)
           data = obj.obj_particles;
           subplot(nplot,1,iplot)
           colo='gbkrm';
           for mat = 1:size(data.Imaterial,1)
                I = data.Imaterial(mat,1):data.Imaterial(mat,2);
                plot(data.Xj(I),y(I),'o','color',colo(mod(mat,4)+1),'MarkerFaceColor','auto'); 
                hold on;
           end
           hold off;
           title([title_name,', t=',num2str(time),' N= ',num2str(data.N)])
           xlim([data.Omega(1) data.Omega(2)]);
           iplot = iplot+1;            
        end            
        %%
        function plot_data2D(obj,t)
           data = obj.obj_particles;
           %% some flags:
           draw_connectivity = false;
           draw_cells = false;
           mark_point = false;
            
           %% draw points
           colo='gbkrm';
           % each index-set(material) with a seperate color
           for mat = 1:size(data.Imaterial,1)
                I = data.Imaterial(mat,1):data.Imaterial(mat,2);
                plot(data.Xj(I,1),data.Xj(I,2),'o','color',colo(mod(mat,4)+1))  %all particles
                hold on;
           end
           axis equal
           rectangle('Position',[data.Omega(1,1),data.Omega(2,1),data.Omega(1,2),data.Omega(2,2)]); %boundary
           title(['position; t=',num2str(t),' N= ',num2str(data.N)])

           %% draw connectivity
           if draw_connectivity               
               colo='gbkrm';
               AiX=1:10:data.N;
               for kk = 1:length(AiX)
                   iX= AiX(kk);
                    neigh_iX = [data.pij(data.pij(:,1)==iX,2);data.pij(data.pij(:,2)==iX,1)]; %neighbours of iX
                    for k=1:size(neigh_iX,1)
                       plot([data.Xj(iX,1);data.Xj(neigh_iX(k),1)],...
                            [data.Xj(iX,2);data.Xj(neigh_iX(k),2)],colo(mod(kk,4)+1));
                    end                    
               end
           end
           
           %% draw cells
           if draw_cells
               %horizontal lines
               for k=1:data.Nc(2)
                  y = data.Omega(2,1) + (k/(data.Nc(2)))*(data.Omega(2,2)-data.Omega(2,1));
                  plot([data.Omega(1,1),data.Omega(1,2)],[y,y],':k')
               end               
               %vertical lines
               for k=1:data.Nc(1)
                  x = data.Omega(1,1) + (k/(data.Nc(1)))*(data.Omega(1,2)-data.Omega(1,1));
                  plot([x,x],[data.Omega(2,1),data.Omega(2,2)],'-k')
               end
           end           
           %% mark one particle
           if mark_point           
               iX=1;%Iin(end)+floor(size(Iboun,1)/2)+1;
               plot(data.Xj(iX,1),data.Xj(iX,2),'rx');  
               % draw cutoff radius
               % r=eta2*data.h;
               % rectangle('Position',[data.Xj(iX,1)-r,data.Xj(iX,2)-r,2*r,2*r],'Curvature',[1,1])
               piX = [find(data.A_pij(:,1)==iX);find(data.A_pij(:,2)==iX)]; %edges of iX
              % niX = [A_pij(A_pij(:,1)==iX,2);A_pij(A_pij(:,2)==iX,1)]; %neighbours of iX

              %todo: update textbox instead of drawing a new one!!!
               annotation('textbox',...
                    [0.15 0.65 0.8 0.25],...
                    'String',{['min(rij) = ' num2str(min(data.A_rij(piX,:))),' (all: ',num2str(min(data.A_rij(:,:))),')'],...
                              ['max(rij) = ' num2str(max(data.A_rij(piX,:))),' (all: ',num2str(max(data.A_rij(:,:))),')'],...
                              ['pj = ' num2str(data.pj(iX)), ' min: ',...
                                    num2str(min(data.pj)),' max: ',num2str(max(data.pj))],...
                              ['Fint = ' num2str(data.F_int(iX,:))],...
                              ['Fdiss = ' num2str(data.F_diss(iX,:))],...
                              ['F ST = ' num2str(data.F_ST(iX,:))]},...
                    'FontSize',10,...
                    'FontName','Arial',...
                    'LineStyle','--',...
                    'EdgeColor',[1 1 0],...
                    'LineWidth',2,...
                    'BackgroundColor',[0.9  0.9 0.9],...
                    'Color',[0.84 0.16 0]);
           end
           
           hold off;
           drawnow
        end
        %%
        function plot_trisurf(obj,time)  
            data = obj.obj_particles;
            x=data.Xj(:,1);
            y=data.Xj(:,2);
            xy = complex(x,y);
            t=delaunay(x,y);
            z=data.pj;
            e = abs(xy(t) - xy(circshift(t,-1,2)));
            goodt = all(e < 2*data.h,2); %dxmedian(e(:)),2);
            t2 = t(goodt,:);
            trisurf(t2,x,y,z)
            %trimesh(t2,x,y,z)

            shading interp
            caxis([-3 3])
            colorbar 
            colormap jet
           % shading flat
            view(2)
          %  axis('equal')
       %     hold off;
            hold on
            rectangle('Position',[0,0,data.Omega(1),data.Omega(2)]); %boundary
            title(['t=',num2str(time),' N= ',num2str(data.N)])
            hold off;
       
            drawnow;
        end
        %%
        function plot_patches(obj,time)
            data = obj.obj_particles;
            x=data.Xj(:,1);
            y=data.Xj(:,2);
            a=data.h;
            z=data.pj;
            opacity = 0.3;
            clf
            transparentScatter(obj,x,y,z,a,opacity);
            %scatter(x,y,a,z,'filled'); % no alpha possible
             caxis([-1 1])
            colorbar 
            colormap jet
            
            hold on
            rectangle('Position',[data.Omega(1,1),data.Omega(2,1),data.Omega(1,2),data.Omega(2,2)]); %boundary
            title(['t=',num2str(time),' N= ',num2str(data.N)])
            axis('equal')

            hold off;
                   
            drawnow;
        end
        %%
        function transparentScatter(~, x,y,z,sizeOfCirlce,opacity)
                t= 0:pi/10:2*pi;
                rep_x = repmat(x',[size(t,2),1]);
                rep_y = repmat(y',[size(t,2),1]);
                rep_z = repmat(z',[size(t,2),1]);
                rep_t = repmat(t',[ 1, size(x,1)]);

                scatterPoints = patch((sizeOfCirlce*sin(rep_t)+ rep_x),...
                                      (sizeOfCirlce*cos(rep_t)+rep_y),...
                                       rep_z,'edgecolor','none');
                alpha(scatterPoints,opacity);

        end
        %%

        %% %%% In/out %%% %%
        function read(obj,filename)
            disp(['read ',filename]);

            filename='in-hvi.h5';
            time = 0;
            mat2h5(obj,filename,time)
            error('not finished');
        end
        
        %ToDo:
        function mat2h5(~,filename,time)
  
            function x=readVariable(x_name,filename,time)
                 x=h5read(filename,['/',time,'/',x_name]);            
            end
            time_str = num2str(time);
             

            Gamma = readVariable('Gamma',filename,time_str);
            Gmod = readVariable('Gmod',filename,time_str);
            S = readVariable('S',filename,time_str);
            Y0 = readVariable('Y0',filename,time_str);
            c = readVariable('c',filename,time_str);
            c0 = readVariable('c0',filename,time_str);
            e = readVariable('e',filename,time_str);
            m = readVariable('m',filename,time_str);
            p = readVariable('p',filename,time_str);
            phi = readVariable('phi',filename,time_str);
            rho = readVariable('rho',filename,time_str);
            rho0 = readVariable('rho0',filename,time_str);
            tauXX = readVariable('tauXX',filename,time_str);
            tauXY = readVariable('tauXY',filename,time_str);
            tauYY = readVariable('tauYY',filename,time_str);
            u = readVariable('u',filename,time_str);
            v = readVariable('v',filename,time_str);
            x = readVariable('x',filename,time_str);
            y = readVariable('y',filename,time_str);

            keyboard

            %% in
            % Gamma, Gmod, S, Y0, c, c0, e, m, p, phi, rho, rho0, 
            % tauXX, tauXY, tauYY,
            % u, v, x,y

            %% out
            % Gamma, Gmod, O, S, Vol, Y0, c, c0, e, eh, epsdotXX, epsdotXY,
            % et, h , m, nn, p ,phi, rho, rho0, rhoh, rhot, rotdotXY,
            % stationary, t, tauXX, tauXXpl, tauXXrho2, tauXY, tauXYpl,
            % tauXYrho2, tauYY, tauYYpl,tauYYrho2, taudotXX, taudotXY, 
            % taudotYY, tmp, u, uh, ups, ut, v, vh, vt, x ,y , 


        end
        

    end

    
end

