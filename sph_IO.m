classdef sph_IO < handle
    %IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %input
        
        %output       
        output_name
        write_data
        
        %plot
        mfigure
        plotstyle
        plot_dt
        
        fixaxes
        
        %movie
        save_as_movie
        movie_name
        vidObj
       
        %conservation variables
        con_dt
        con_mass
        con_momentum
        con_energy
        
    end
    
    methods
        
        %% %%% constructor
        function obj = sph_IO(obj_scen)
              % some IO properties
              obj.save_as_movie = obj_scen.save_as_movie;
              obj.movie_name    = get_movie_name(obj_scen);
              obj.plot_dt       = obj_scen.plot_dt;
              obj.plotstyle     = obj_scen.plotstyle;
                                         
              %read data from file and save in obj_scen
              if obj_scen.read_data
                 obj.read_hdf5(obj_scen);
              end    
              
              %output
              obj.write_data  = obj_scen.write_data;
              obj.output_name = obj_scen.output_name;
              if obj.write_data %create file
                  hdf5write(obj.output_name, '/name', 0); % ToDo: better solution?
              end   
              
              %conservation variables              
              obj.con_mass = [];
              obj.con_momentum =[];
              obj.con_dt = [];
              obj.con_energy =[];
              obj.fixaxes = obj_scen.fixaxes;
        end
        
        %% %%% general functions
        function initialize(obj)
            %% output

            %% plot
            if ~strcmp(obj.plotstyle,'')  %plot only, if plotstlye is specified
                  obj.mfigure = figure;
            end

            
            %% movie
            if obj.save_as_movie
                set(gca,'DataAspectRatio',[1,1,1]);
                obj.vidObj     = VideoWriter(obj.movie_name);
                open(obj.vidObj);
            end            
            
        end
        %%
        function do (obj, obj_particles)
            %% plotting
            if (mod(obj_particles.t,obj.plot_dt) < obj_particles.dt) %ToDo: good for variable timestepping?
                plot_data (obj,obj_particles);
                if obj.save_as_movie
                    currFrame = getframe(obj.mfigure);
                    writeVideo(obj.vidObj,currFrame);
                end
               
                %write output
                if obj.write_data
                     write_hdf5(obj,obj_particles)
                end
            end  
            %check for conservation of mass and momentum
            checkConservation (obj,obj_particles);
 

        end       
        %%
        function finalize (obj)
            %% close output
            
            %% movie
            if obj.save_as_movie
              close(obj.vidObj);
              disp (['movie saved as ',obj.movie_name]);
            end   
            plot_conservation = true;
            mom_norm = sum(obj.con_momentum.*obj.con_momentum,2).^0.5;
            if plot_conservation
                %% plot conservation variables
                fig=figure;
                subplot(3,1,1)
                plot(obj.con_dt,obj.con_mass)
                title('evolution of mass');
                xlabel('t'); ylabel('mass')
                subplot(3,1,2)            
                % norm(momentum)
                plot(obj.con_dt,mom_norm);            
                hold on;
                %components (in 2d)
                dim=size(obj.con_momentum,2);            
                if dim > 1
                    for i = 1: dim
                        plot(obj.con_dt,obj.con_momentum(:,i)) %x
                    end
                    legend('norm','x','y');
                end
                title('evolution of momentum');
                xlabel('t'); ylabel('momentum ')
                hold off;
                subplot(3,1,3)
                plot(obj.con_dt,obj.con_energy)
                %move figure to the left side
                figpos=fig.Position;
                figpos(1)=0;
                fig.Position =figpos;
            end
            disp (['## relative dissipation of momentum = ',...
                num2str(abs((mom_norm(end)-mom_norm(1))/mom_norm(1))),' ##']);
        end
        %%
        function checkConservation (obj,obj_particles)
            data=obj_particles;
            %time
            obj.con_dt = [obj.con_dt;
                          data.t];
            
            % conservation of mass
            massj = data.rhoj.*data.Vj;
            obj.con_mass = [obj.con_mass;...
                            sum(massj)];
            % conservation of momentum
            momj  = massj*ones(1,data.dim) .* data.vj;
            obj.con_momentum = [obj.con_momentum;...
                                sum(momj,1)];
                            
            % energy after Modave 1/2mv^2 + pV
%             ej = 0.5 * massj .* sum(data.vj.^2,2) + data.Vj.*data.pj.^2;
%             obj.con_energy =[obj.con_energy;
%                              sum(ej)];
            % change of energy
            ej = data.vj(data.Iin)'*data.F_total(data.Iin) + ...
                 sum(data.pj(data.Iin)./data.rhoj(data.Iin).^2 ...
               .* data.drhoj(data.Iin) .*data.mj(data.Iin));
            obj.con_energy =[obj.con_energy;
                              sum(ej)];
        end

        
        %% %%%  Plotting functions %%% %%
        function plot_data(obj,obj_particles)
           if ~strcmp(obj.plotstyle,'')  %plot only, if plotstlye is specified
               figure(obj.mfigure);
               if obj_particles.dim == 1               
                   plot_data1D(obj,obj_particles);
               elseif obj_particles.dim == 2
                   if strcmp(obj.plotstyle,'scatter')
                        plot_scatter2D(obj,obj_particles);
                   elseif strcmp (obj.plotstyle,'trisurf')
                        plot_trisurf(obj,obj_particles)
                   elseif strcmp (obj.plotstyle,'patches')
                        plot_patches(obj,obj_particles)
                   else
                       error([obj.plotstyle, '- plotstyle is not supported']);
                   end
               end 
           end
        end          
        %%
        function plot_data1D(obj,obj_particles)
            %
            function draw_nr_bc_area(data,y)
                %%-- draw non-reflecting boundary condition
               if ~isempty(data.damping_area)
                   a = axis;
                   b = data.omega_nr_bc;
                   rectangle('Position',[b(1), a(3), b(2)-b(1), a(4)-a(3)],...
                             'EdgeColor',[0,1,0.5]);
                   axis(a)
                   kb = logical((data.Xj>b(1)) .* (data.Xj<b(2))); %particels in boundary layer
                   % mark particles
                   plot(data.Xj(kb),y(kb),'xr');
               end
               %%--                               
            end
            %
            function iplot = create_a_supplot_1d(data,y,title_name,nplot,iplot,yaxes)
                   subplot(nplot,1,iplot)
                   colo='gbkrm';
                   time = data.t;
                   %each material gets his own color
                   for mat = 1:size(data.Imaterial_with_boun,1)
                        I = data.Imaterial_with_boun(mat,1):data.Imaterial_with_boun(mat,2);
                        plot(data.Xj(I),y(I),'o','color',colo(mod(mat,4)+1),'MarkerFaceColor','auto'); 
                        hold on;
                   end                     
                   xlim([data.Omega(1) data.Omega(2)]);
                   if ~isempty(yaxes)
                        ylim(yaxes);
                   end
                   draw_nr_bc_area(data,y) %draw non-reflecting boundary condition
                   hold off;
                   title([title_name,', t=',num2str(time),' N= ',num2str(data.N)])
                   iplot = iplot+1;                               
            end 
            
            data  = obj_particles;
            nplot = length(obj.plotstyle);
            iplot = 1;
            if ~isempty(strfind(obj.plotstyle,'x'));
                create_a_supplot_1d(data,zeros(data.N,1),'position',nplot,iplot,obj.fixaxes.x);
                iplot=iplot+1;     
            end
            if ~isempty(strfind(obj.plotstyle,'v'));
                iplot = create_a_supplot_1d(data,data.vj,'velocity',nplot,iplot,obj.fixaxes.v);
            end
            if ~isempty(strfind(obj.plotstyle,'p'));
                iplot = create_a_supplot_1d(data,data.pj,'pressure',nplot,iplot,obj.fixaxes.p);
            end
            if ~isempty(strfind(obj.plotstyle,'d'));
                iplot = create_a_supplot_1d(data,data.rhoj,'density',nplot,iplot,obj.fixaxes.d);
            end
            if ~isempty(strfind(obj.plotstyle,'f'));     
                subplot(nplot,1,iplot)
                bar(data.Xj,[data.F_int,data.F_diss,data.F_diss_art,data.F_ST]); title('forces')
                xlim([data.Omega(1) data.Omega(2)]);
                legend('F_{int}','F_{diss}','F_{diss-art}','F_{ST}');
            end
            drawnow
        end
        %%
        function plot_scatter2D(~,obj_particles)
           data = obj_particles;
           t = obj_particles.t;
           %% some flags:
           draw_connectivity = false;
           draw_cells = false;
           mark_point = false;
            
           %% draw points
           colo='gbkrm';
           % each index-set(material) with a seperate color
           for mat = 1:size(data.Imaterial_with_boun,1)
                I = data.Imaterial_with_boun(mat,1):data.Imaterial_with_boun(mat,2);
                plot(data.Xj(I,1),data.Xj(I,2),'o','color',colo(mod(mat,4)+1))  %all particles
                hold on;
           end
           axis equal
           rectangle('Position',[data.Omega(1,1),data.Omega(2,1),...
               data.Omega(1,2)-data.Omega(1,1),data.Omega(2,2)-data.Omega(2,1)]); %boundary
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
                              ['F_{int} = ' num2str(data.F_int(iX,:))],...
                              ['F_{diss} = ' num2str(data.F_diss(iX,:))],...
                              ['F_{ST} = ' num2str(data.F_ST(iX,:))]},...
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
        function plot_trisurf(~,obj_particles)  
            data = obj_particles;
            time = data.t;
            x    = data.Xj(:,1);
            y    = data.Xj(:,2);
            xy   = complex(x,y);
            t    = delaunay(x,y);
            z    = data.pj;
            e    = abs(xy(t) - xy(circshift(t,-1,2)));
            goodt = all(e < 2*max(data.h),2); %dxmedian(e(:)),2);
            t2   = t(goodt,:);
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
            rectangle('Position',[data.Omega(1,1),data.Omega(2,1),...
                     data.Omega(1,2)-data.Omega(1,1),data.Omega(2,2)-data.Omega(2,1)]); %boundary
            title(['t=',num2str(time),' N= ',num2str(data.N)])
            hold off;
       
            drawnow;
        end
        %%
        function plot_patches(obj,obj_particles)
            
            function transparentScatter(x,y,z,sizeOfCirlce,opacity)
                t= 0:pi/10:2*pi;
                rep_x = repmat(x',[size(t,2),1]);
                rep_y = repmat(y',[size(t,2),1]);
                rep_r = repmat(sizeOfCirlce',[size(t,2),1]);
                rep_z = repmat(z',[size(t,2),1]);
                rep_t = repmat(t',[ 1, size(x,1)]);

                scatterPoints = patch((rep_r.*sin(rep_t)+ rep_x),...
                                      (rep_r.*cos(rep_t)+rep_y),...
                                       rep_z,'edgecolor','none');
                alpha(scatterPoints,opacity);
             end
            data = obj_particles;
            time = data.t;
            x=data.Xj(:,1);
            y=data.Xj(:,2);
            a=data.hj;
            z=data.pj;  %visualize the pressure
            opacity = 0.3;
            clf
            transparentScatter(x,y,z,a,opacity);
            if ~isempty(obj.fixaxes.p)
                caxis(obj.fixaxes.p)
            end
            colorbar 
            colormap jet
            
            hold on
            rectangle('Position',[data.Omega(1,1),data.Omega(2,1),...
               data.Omega(1,2)-data.Omega(1,1),data.Omega(2,2)-data.Omega(2,1)]); %boundary
                title(['t=',num2str(time),' N= ',num2str(data.N)])
            axis('equal')

            hold off;
                   
            drawnow;
        end
        

        %% %%% In/out %%% %%
        %%
        function read_hdf5(~,obj_scen)
  
            function x=readVariable(x_name,filename,time)
                 x=h5read(filename,['/',time,'/',x_name]);            
            end
            time = 0;
            time_str = num2str(time);
            filename =obj_scen.input_name;
            %position
            obj_scen.Xj = [readVariable('x',filename,time_str),...
                           readVariable('y',filename,time_str)];
            %velocity                       
            obj_scen.vj = [readVariable('u',filename,time_str),...
                           readVariable('v',filename,time_str)];
    
            %speed of sound
            %c = readVariable('c',filename,time_str);
            obj_scen.c0j = readVariable ('c0',filename,time_str);
            obj_scen.cj = readVariable('c',filename,time_str);

            %density            
            %rho = readVariable('rho',filename,time_str);
            obj_scen.rho0j = readVariable('rho0',filename,time_str);
            obj_scen.rhoj = readVariable('rho',filename,time_str);
     
            %mass
            obj_scen.mj = readVariable('m',filename,time_str);
            
            %ToDo: improve this
            obj_scen.dx  = abs(obj_scen.Xj(1,1)-obj_scen.Xj(2,1));
            N = size(obj_scen.Xj,1);
            obj_scen.Iin = (1:N)';
            obj_scen.Imaterial = [1,N]; 
          %  Gamma = readVariable('Gamma',filename,time_str);
          %  Gmod = readVariable('Gmod',filename,time_str);
          %  S = readVariable('S',filename,time_str);
          %  Y0 = readVariable('Y0',filename,time_str);
          %  e = readVariable('e',filename,time_str);
          %  p = readVariable('p',filename,time_str);
          %  phi = readVariable('phi',filename,time_str);
          %  tauXX = readVariable('tauXX',filename,time_str);
          %  tauXY = readVariable('tauXY',filename,time_str);
          %  tauYY = readVariable('tauYY',filename,time_str);

        end
        %%
        function write_hdf5(obj,obj_particle)
            filename = obj.output_name;
            time = obj_particles.t;
            function writeVariable(filename,time,name,data)
                group = ['/',num2str(time),'/',name];
                hdf5write(filename, group, data, 'WriteMode', 'append');
            end
            writeVariable(filename,time,'x',obj_particle.Xj(:,1));
            writeVariable(filename,time,'y',obj_particle.Xj(:,2));
            writeVariable(filename,time,'u',obj_particle.vj(:,1));
            writeVariable(filename,time,'v',obj_particle.vj(:,2));
            writeVariable(filename,time,'c',obj_particle.cj);
            writeVariable(filename,time,'c0',obj_particle.cj0);
            writeVariable(filename,time,'rho',obj_particle.rhoj);
            writeVariable(filename,time,'rho0',obj_particle.rho0j);
            writeVariable(filename,time,'m',obj_particle.mj);                        
            
            %% out
            % Gamma, Gmod, O, S, Vol, Y0, c, c0, e, eh, epsdotXX, epsdotXY,
            % et, h , m, nn, p ,phi, rho, rho0, rhoh, rhot, rotdotXY,
            % stationary, t, tauXX, tauXXpl, tauXXrho2, tauXY, tauXYpl,
            % tauXYrho2, tauYY, tauYYpl,tauYYrho2, taudotXX, taudotXY, 
            % taudotYY, tmp, u, uh, ups, ut, v, vh, vt, x ,y , 
        end
    end

    
end

