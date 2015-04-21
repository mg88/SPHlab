classdef sph_IO < handle
    %IO classes for the sph-simulation tool
    
    properties
        %input
        
        %output       
        output_name
        write_data
        
        %plot
        mfigure       %figurehandle
        plot_style
        plot_param
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
              obj.plot_style    = obj_scen.plot_style;
              obj.plot_param    = obj_scen.plot_param;
                                         
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
            if ~strcmp(obj.plot_param,'')  %plot only, if plotstlye is specified
                 obj.mfigure = figure;
                 set(gca,'DataAspectRatio',[1,1,1]);
                 if length(obj.plot_param) > 3
                        obj.mfigure.Units = 'normalized';
                        figpos=obj.mfigure.Position;
                        figpos(1)=0.1;
                        figpos(3)=0.8;
                        figpos(2)=0.1;
                        figpos(4)=0.8;
                        obj.mfigure.Position =figpos;
                   elseif length(obj.plot_param) > 1 
                        obj.mfigure.Units = 'normalized';
                        figpos=obj.mfigure.Position;
                        figpos(1)=0.1;
                        figpos(3)=0.8;
                        obj.mfigure.Position =figpos;
                 end
            end

            
            %% movie
            if obj.save_as_movie
%                 set(gca,'DataAspectRatio',[1,1,1]);
%                 obj.mfigure.Units = 'normalized';
%                 obj.mfigure.Position = [0 0 1 1];
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
                title('change of energy');
                %move figure to the left side
                figpos=fig.Position;
                figpos(2)=0;
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
            massj = data.rhoj(data.Iin).*data.Vj(data.Iin);
            obj.con_mass = [obj.con_mass;...
                            sum(massj)];
            % conservation of momentum
            momj  = massj*ones(1,data.dim) .* data.vj(data.Iin,:);
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
            
            function plot_scatter(x,dat,mat)
               colo='gbkrm';
               %each material gets his own color
               for m = 1:size(mat,1)
                    I = mat(m,1):mat(m,2);
                    plot(x(I),dat(I),'o','color',colo(mod(m,4)+1),'MarkerFaceColor','auto'); 
               end                     
            end    
             
            function plot_force(data)
                bar(data.Xj(data.Iin),...
                    [data.F_int(data.Iin),data.F_diss(data.Iin),data.F_diss_art(data.Iin),data.F_ST(data.Iin)]);
                title('forces')
                xlim([data.Omega(1) data.Omega(2)]);
                legend('F_{int}','F_{diss}','F_{diss-art}','F_{ST}');
            end
           
            function plot_patches(x,y,z,sizeOfCirlce,opacity)
                tt= 0:pi/10:2*pi;
                rep_x = repmat(x',[size(tt,2),1]);
                rep_y = repmat(y',[size(tt,2),1]);
                rep_r = repmat(sizeOfCirlce',[size(tt,2),1]);
                rep_z = repmat(z',[size(tt,2),1]);
                rep_tt = repmat(tt',[ 1, size(x,1)]);
                scatterPoints = patch((rep_r.*sin(rep_tt)+ rep_x),...
                                      (rep_r.*cos(rep_tt)+rep_y),...
                                       rep_z,'edgecolor','none');
                alpha(scatterPoints,opacity);
            end 

            function plot_trisurf(x,y,z,max_r)  
                xy   = complex(x,y);
                tri    = delaunay(x,y);
                e    = abs(xy(tri) - xy(circshift(tri,-1,2)));
                goodt = all(e < max_r,2); %dxmedian(e(:)),2);
                t2   = tri(goodt,:);
                trisurf(t2,x,y,z)
                %trimesh(t2,x,y,z)
            end
                     
            function dat_max=plot_field(x,dat)
               quiver(x(:,1),...
                         x(:,2),...
                         dat(:,1),...
                         dat(:,2));
               dat_norm = sum(dat.^2,2).^0.5;
               dat_max = max(dat_norm);
            end
            
           figure(obj.mfigure);
           nplot = length(obj.plot_param);
           if nplot <=3    
               nxplot = nplot;
               nyplot =1;
           else
               nxplot = ceil(nplot/2);
               nyplot = 2;
           end
           iplot=1;
           x   = obj_particles.Xj;
           t   = obj_particles.t;
           mat = obj_particles.Imaterial_with_boun;
           title_additive = ['; t=',num2str(t),'; N= ',num2str(obj_particles.N)];
           clf %clear figure
           for para = obj.plot_param;
               subplot(nyplot,nxplot,iplot);
               hold on;
               if para == 'x'
                   dat = zeros(obj_particles.N,1);
                   name = 'position';
                   limaxes = obj.fixaxes.x;
                   style   = obj.plot_style.x;
               elseif para == 'p'
                   dat = obj_particles.pj;
                   name = 'pressure';
                   limaxes = obj.fixaxes.p;
                   style   = obj.plot_style.p;
               elseif para == 'd'
                   name = 'density';
                   dat = obj_particles.rhoj;
                   limaxes = obj.fixaxes.d;                   
                   style   = obj.plot_style.d;
               elseif para == 'v'
                   name = 'velocity';
                   dat = obj_particles.vj;
                   limaxes = obj.fixaxes.v;
                   style   = obj.plot_style.v;
               elseif para == 'f'
                   name = 'F-total'; %todo plot components
                   dat = obj_particles.F_total;
                   limaxes = obj.fixaxes.f;
                   style   = obj.plot_style.f;
               elseif para == 'e'
                   name = 'energy';
                   dat = obj_particels.ej;
                   limaxes = obj.fixaxes.e;
                   style   = obj.plot_style.e;
               else
                   warning('variable is not implemented for plotting');
                   keyboard
               end

               %plot:
               if obj_particles.dim == 1   
                   if para == 'f'
                       plot_force(obj_particles);
                   else
                       plot_scatter(x,dat,mat); 
                       
                       %mark mirror particle
                       plot(obj_particles.Xj(obj_particles.bc.mirrorParticlesj),dat(obj_particles.bc.mirrorParticlesj),'xr');
                       
                   end
                   if ~isempty(limaxes)
                     ylim(limaxes);
                   end                   
               elseif obj_particles.dim == 2
                   if para == 'x'
                      plot_scatter(x(:,1),x(:,2),mat); 
                      axis equal
                   elseif any(para == 'pde') %perssure, density, energy -> scalar
                        if strcmp (style,'trisurf')
                            plot_trisurf(x(:,1),x(:,2),dat,2*max(obj_particles.hj));
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            view([-1,-1,1]);
                            %view(2)
                            shading interp
                            if ~isempty(limaxes)
                                 zlim(limaxes)
                            end
                            colorbar 
                            colormap jet
                        elseif strcmp (style,'patches')                            	
                            opacity = 0.3;
                            sizeOfCirlce = obj_particles.hj;
                            plot_patches(x(:,1),x(:,2),dat,sizeOfCirlce,opacity)                 
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            colorbar 
                            colormap jet
                            axis equal
                        elseif strcmp (style,'plot3')
                            plot3(x(:,1),x(:,2),dat,'o');
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                                 zlim(limaxes)
                            end
                            view([0,-1,0]);
                        else
                            error([style, '- plotstyle is not supported']);
                        end
                   elseif any(para == 'vf') % velocity, forces -> field
                           dat_max = plot_field(x(obj_particles.Iin,:),dat(obj_particles.Iin,:));
                           axis equal
                           title_additive= [title_additive,'; max|',para,'|=',num2str(dat_max)];
                   else
                       error([obj.plotstyle, '- plotstyle is not supported']);
                   end                   
               end 
               title([name,title_additive]);
               title_additive='';
               %% plot additional info
               plot_geometry(obj,obj_particles);
               iplot = iplot+1;
               hold off;
           end
           drawnow
        end          
        %%
        function plot_geometry(~,obj_p)
           mark_point = false;
           draw_cells = false;
           draw_connectivity = false;                %only 2d

           if obj_p.dim == 1
               % Omega
               xlim( [obj_p.Omega(1), obj_p.Omega(2)]);
               
                %%-- draw damping boundary condition
               for boun = obj_p.bc
                   if ~isempty(boun.damping_area)
                       a = axis;
                       b = boun.damping_area;
                       rectangle('Position',[b(1), a(3), b(2)-b(1), a(4)-a(3)],...
                                 'EdgeColor',[0,1,0.5]);
                       axis(a);
                   end
               end
                
           elseif obj_p.dim == 2
               
               %axis('equal')                   
               %Omega:
               rectangle('Position',[obj_p.Omega(1,1),obj_p.Omega(2,1),...
                   obj_p.Omega(1,2)-obj_p.Omega(1,1),obj_p.Omega(2,2)-obj_p.Omega(2,1)]); %boundary
               hold on


               %% draw connectivity
               if draw_connectivity               
                   colo='gbkrm';
                   AiX=1:500:obj_p.N;
                   for kk = 1:length(AiX)
                       iX= AiX(kk);
                        neigh_iX = [obj_p.pij(obj_p.pij(:,1)==iX,2);obj_p.pij(obj_p.pij(:,2)==iX,1)]; %neighbours of iX
                        for k=1:size(neigh_iX,1)
                           plot([obj_p.Xj(iX,1);obj_p.Xj(neigh_iX(k),1)],...
                                [obj_p.Xj(iX,2);obj_p.Xj(neigh_iX(k),2)],colo(mod(kk,4)+1));
                        end                    
                   end
               end



           end
           
           %% draw cells
           if draw_cells
               disp 'draw cell-structure';
               if obj_p.dim == 2
                   %horizontal lines
                   for k=1:obj_p.Nc(2)
                      y = obj_p.Omega(2,1) + (k/(obj_p.Nc(2)))*(obj_p.Omega(2,2)-obj_p.Omega(2,1));
                      plot([obj_p.Omega(1,1),obj_p.Omega(1,2)],[y,y],':k')
                   end 
                   vline = [obj_p.Omega(2,1),obj_p.Omega(2,2)]
               else
                   a=axis;
                   vline =[a(3) ,a(4)];
               end
               %vertical lines
               for k=1:obj_p.Nc(1)
                  x = obj_p.Omega(1,1) + (k/(obj_p.Nc(1)))*(obj_p.Omega(1,2)-obj_p.Omega(1,1));
                  plot([x,x],vline,'-k')
               end
           end  
           
                      %% mark one particle
           if mark_point           
               iX=1;%Iin(end)+floor(size(Iboun,1)/2)+1;
               if obj_p.dim == 1
                   plot(obj_p.Xj(iX,1),0,'rx');  
               else
                   plot(obj_p.Xj(iX,1),obj_p.Xj(iX,2),'rx');  
               end
               piX = [find(obj_p.pij(:,1)==iX);find(obj_p.pij(:,2)==iX)]; %edges of iX
              % niX = [pij(pij(:,1)==iX,2);A_pij(pij(:,2)==iX,1)]; %neighbours of iX

              %todo: update textbox instead of drawing a new one!!!
               annotation('textbox',...
                    [0.15 0.65 0.8 0.25],...
                    'String',{['min(rij) = ' num2str(min(obj_p.rij(piX,:))),' (all: ',num2str(min(obj_p.rij(:,:))),')'],...
                              ['max(rij) = ' num2str(max(obj_p.rij(piX,:))),' (all: ',num2str(max(obj_p.rij(:,:))),')'],...
                              ['pj = ' num2str(obj_p.pj(iX)), ' min: ',...
                                    num2str(min(obj_p.pj)),' max: ',num2str(max(obj_p.pj))],...
                              ['F_{int} = ' num2str(obj_p.F_int(iX,:))],...
                              ['F_{diss} = ' num2str(obj_p.F_diss(iX,:))],...
                              ['F_{ST} = ' num2str(obj_p.F_ST(iX,:))]},...
                    'FontSize',10,...
                    'FontName','Arial',...
                    'LineStyle','--',...
                    'EdgeColor',[1 1 0],...
                    'LineWidth',2,...
                    'BackgroundColor',[0.9  0.9 0.9],...
                    'Color',[0.84 0.16 0]);
           end
           
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

