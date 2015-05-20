classdef sph_IO < handle
    %IO classes for the sph-simulation tool
    
    properties
        % evaluation points
        x_eval
        dx_eval
        Nss  %super sampling
        %input
        
        %output       
        output_name
        write_data
        save_dt
        
        %plot
        mfigure       %figurehandle
        plot_style
        plot_quantity
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
        con_momentum_diss
        con_e_kin
        con_e_kin_diss
        con_e_pot
        con_e_pot_diss
        
        %some internal flags/variables
        t_last_plot
        t_last_save
     
        %exact solution
        exact_sol
    end
    
    methods
        
        %% %%% constructor
        function obj = sph_IO(obj_scen)
            if nargin == 1 %else, use only the functions
              %% some IO properties
              obj.save_as_movie = obj_scen.save_as_movie;
              obj.movie_name    = get_movie_name(obj_scen);
              obj.save_dt       = obj_scen.save_dt;
              obj.plot_dt       = obj_scen.plot_dt;
              obj.plot_style    = obj_scen.plot_style;
              obj.plot_quantity = obj_scen.plot_quantity;                                          
              obj.exact_sol     = obj_scen.exact_sol;                                          

              %% gridpoints to evaluate
              if obj_scen.Neval > 0
                  [obj.x_eval,obj.dx_eval] = obj_scen.rectangle(obj_scen.Omega,obj_scen.Neval);
              else
                  obj.x_eval  = [];
                  obj.dx_eval = [];              
              end
              obj.Nss     = obj_scen.Nss;

              
              %% output
              obj.write_data  = obj_scen.write_data;
              obj.output_name = obj_scen.output_name;
              if obj.write_data %create file
                  hdf5write([obj.output_name,'.h5'], '/name', 0); % ToDo: better solution?
                  save([obj.output_name,'_scen.mat'], 'obj_scen');
              end   
              
              %% conservation variables              
              obj.con_mass = [];
              obj.con_momentum =[];
              obj.con_momentum_diss =[];
              obj.con_dt = [];
              obj.con_e_kin =[];
              obj.con_e_kin_diss =[];
              obj.con_e_pot =[];
              obj.con_e_pot_diss =[];
              obj.fixaxes = obj_scen.fixaxes;
              
              %%
              obj.t_last_plot = -inf;
              obj.t_last_save = -inf;
            end
        end
        
        %% %%% general functions
        function initialize(obj)

            %% plot
            if ~strcmp(obj.plot_quantity,'')  %plot only, if plotstlye is specified
                 obj.mfigure = figure;
                 set(gca,'DataAspectRatio',[1,1,1]);
                 if strcmp(get(0,'DefaultFigureWindowStyle'),'normal') %only when not docked
                     if length(obj.plot_quantity) > 3
                            obj.mfigure.Units = 'normalized';
                            figpos=obj.mfigure.Position;
                            figpos(1)=0.1;
                            figpos(3)=0.8;
                            figpos(2)=0.1;
                            figpos(4)=0.8;
                            obj.mfigure.Position =figpos;
                       elseif length(obj.plot_quantity) > 1 
                            obj.mfigure.Units = 'normalized';
                            figpos=obj.mfigure.Position;
                            figpos(1)=0.1;
                            figpos(3)=0.8;
                            obj.mfigure.Position =figpos;
                     end
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
            
            %% output
            
        end
        %%
        function do (obj, obj_particles,force)
            if nargin <3
                force = false;
            end
            %% plotting
            if force||(((obj_particles.t-obj.t_last_plot) > obj.plot_dt) ...
                    && ~isempty(obj.plot_quantity))
                plot_data (obj,obj_particles);
                if obj.save_as_movie
                    currFrame = getframe(obj.mfigure);
                    writeVideo(obj.vidObj,currFrame);
                end         
                obj.t_last_plot = obj_particles.t;
            end  
            
            %% write output
            if (force||((obj_particles.t-obj.t_last_save) > obj.plot_dt)) ...
                            && obj.write_data
                 write_hdf5(obj,obj_particles)
                %name = ['data/out',num2str(obj_particles.dt),'.mat'];               
                % save(name, 'obj_particles');   
                 obj.t_last_save = obj_particles.t;
            end
            
            %check for conservation of mass and momentum
            save_conservation (obj,obj_particles);

        end       
        %%
        function finalize (obj)
            
            %% close output
            
            %% movie
            if obj.save_as_movie
              close(obj.vidObj);
              disp (['movie saved as ',obj.movie_name]);
            end   
            
            %% plot conservation variables
            plot_conservation(obj)

        end
        %%
        function save_conservation (obj,obj_particles)
            data=obj_particles;
            %time
            obj.con_dt = [obj.con_dt;
                          data.t];
            
            % conservation of mass
            massj = data.rhoj(data.Iin).*data.Vj(data.Iin);
            obj.con_mass = [obj.con_mass;...
                            sum(massj)];
                        
            % conservation of momentum (particle inside)
            momj  = massj * ones(1,data.dim) .* (data.vj_half(data.Iin,:));
            obj.con_momentum = [obj.con_momentum;...
                                sum(momj,1)];
            
            % dissipation of momentum  (at half step)          
            if (~isempty(data.Ighost) && data.compGhost) 
                momj_diss_add   = data.IOdata_h.vjghost.*(data.mj(data.Ighost)*ones(1,data.dim));                
                momj_diss_after = data.vj_half(data.Ighost,:).*(data.mj(data.Ighost)*ones(1,data.dim));

                if isempty(obj.con_momentum_diss) %time integration
                    previous = zeros(1,data.dim); 
                else
                    previous = obj.con_momentum_diss(end,:);
                end
                obj.con_momentum_diss =[obj.con_momentum_diss;...
                                        sum(momj_diss_after-momj_diss_add,1)+previous];
            end
            
            %energy:
            if data.dim==2
                norm_vj_sqrt = (data.vj_half(data.Iin,1).^2 + data.vj_half(data.Iin,2).^2);
            else
                norm_vj_sqrt = data.vj_half(data.Iin,1).^2;
            end
            e_kin = 0.5 .* data.mj(data.Iin).* norm_vj_sqrt;
            e_pot = data.mj(data.Iin).*data.ej_half(data.Iin);
            obj.con_e_kin =[obj.con_e_kin;
                                  sum(e_kin)];
            obj.con_e_pot =[obj.con_e_pot;
                            sum(e_pot)];
                        
            % energy dissipation
            if (~isempty(data.Ighost) && data.compGhost)
                %dv = data.vj(data.Ighost) - data.IOdata.vjghost;
                if data.dim==2
                    norm_vj_sqrt = (data.vj_half(data.Ighost,1).^2 + data.vj_half(data.Ighost,2).^2);
                    norm_vj_sqrt_prev = (data.IOdata_h.vjghost(:,1).^2 + data.IOdata_h.vjghost(:,2).^2);
                else
                    norm_vj_sqrt = data.vj_half(data.Ighost,1).^2;
                    norm_vj_sqrt_prev = data.IOdata_h.vjghost(:,1).^2;
                end
                
                e_kin_diss      = 0.5 .* data.mj(data.Ighost).* norm_vj_sqrt;
                e_kin_diss_prev = 0.5 .* data.mj(data.Ighost).* norm_vj_sqrt_prev;                
                e_kin_diss_diff =  e_kin_diss - e_kin_diss_prev;
                
                e_pot_diss      = data.ej(data.Ighost).*data.mj(data.Ighost);
                e_pot_diss_prev = data.IOdata_h.ejghost.*data.mj(data.Ighost);
                e_pot_diss_diff = e_pot_diss - e_pot_diss_prev;
                

                if isempty(obj.con_e_kin_diss)
                    previous_kin = 0;
                    previous_pot = 0;
                else                
                    previous_kin = obj.con_e_kin_diss(end);
                    previous_pot = obj.con_e_pot_diss(end);
                end

                obj.con_e_kin_diss =[obj.con_e_kin_diss;
                                      sum(e_kin_diss_diff)+previous_kin];
                obj.con_e_pot_diss =[obj.con_e_pot_diss;
                                     sum(e_pot_diss_diff)+previous_pot];
            end
        end
        %%
        function plot_conservation (obj)
            fig=figure;
            nplot = 3; 
            iplot = 1;            
            if any(diff(obj.con_mass) ~= 0)  % plot only if change of mass occurs              
                subplot(nplot,1,iplot)
                iplot=iplot+1;
                plot(obj.con_dt,obj.con_mass)
                title('evolution of mass');
                xlabel('t'); ylabel('mass')
            else
                nplot = nplot-1;
            end
            
            %momentum
            subplot(nplot,1,iplot)
            grid on;
            hold on;
            iplot=iplot+1;
            %components (in 2d)
            dim = size(obj.con_momentum,2);            
            if dim ==1
                plot(obj.con_dt,obj.con_momentum); 
                leg = 'I         ';
                %dissipation:
                if ~isempty(obj.con_momentum_diss)
                    plot(obj.con_dt,obj.con_momentum_diss,':');   
                    plot(obj.con_dt,obj.con_momentum_diss+obj.con_momentum,'k-.');   
                    leg = [leg;'I_{diss}  ';'I_{tot}   '];
                end
                mom_norm = abs(obj.con_momentum);           
                
            else
                mom_norm = sum(obj.con_momentum.*obj.con_momentum,2).^0.5;
                plot(obj.con_dt,mom_norm);    
                leg = '|I|       ';
                %components:
                for i = 1: dim
                    plot(obj.con_dt,obj.con_momentum(:,i)) %x
                    leg = [leg;'I_{x',num2str(i),'}    '];
                end
                %dissipation:
                mom_norm_diss = sum(obj.con_momentum_diss.*obj.con_momentum_diss,2).^0.5;
                if ~isempty(mom_norm_diss)
                    plot(obj.con_dt,mom_norm_diss,':');   
                    plot(obj.con_dt,mom_norm_diss+mom_norm,'k-.');   
                    leg = [leg;'|I_{diss}|';'|I_{tot}| '];
                end
            end
            

           
            legend(leg);
            title('evolution of momentum (at half step)');
            xlabel('t'); ylabel('momentum ')
            
            %energy
            subplot(nplot,1,iplot)   
            hold on;
            grid on
            if ~isempty(obj.con_e_pot) && ~isempty(obj.con_e_kin) && ~any(obj.con_e_kin==0)                
                plot(obj.con_dt,obj.con_e_pot,...
                     obj.con_dt,obj.con_e_kin,...
                     obj.con_dt,obj.con_e_pot+obj.con_e_kin,'--');  
                 leg = ['e_{pot}  ';'e_{kin}  ';'e_{tot}  '];
                 if ~isempty(obj.con_e_pot_diss)
                     plot(obj.con_dt,obj.con_e_pot_diss,':',...
                        obj.con_dt,obj.con_e_kin_diss,':',...
                        obj.con_dt,obj.con_e_pot_diss+obj.con_e_kin_diss,'--'); 
                     leg = [leg;'e_{pot,d}';'e_{kin,d}';'e_{tot,d}'];
                     plot(obj.con_dt,obj.con_e_pot+obj.con_e_kin+...
                                     obj.con_e_pot_diss+obj.con_e_kin_diss,'k-.');
                     leg = [leg;'sum(e)   '];
                 end
                 legend (leg);
                 title('energy (at half step)');
            else
                warning('cannot plot energy');
            end
            %move figure to the left side (only when not docked)
            if strcmp(get(0,'DefaultFigureWindowStyle'),'normal')
                figpos=fig.Position;
                figpos(2)=0;
                fig.Position =figpos;
            end
            if ~isempty(mom_norm)
                disp (['## relative dissipation of momentum = ',...
                    num2str(abs((mom_norm(end)-mom_norm(1))/mom_norm(1))),' ##']);
            end
        end
        %% super sampling
        function dat_eval = eval_ss (obj,obj_p, data_name)  %ToDo: also for 2d!         
            %define evaluation points for supersampling
           Nsupersampling = obj.Nss;
           Neval   = size(obj.x_eval,1);
           
           Noff_right = floor(Nsupersampling/2);
           if mod(Nsupersampling,2)==0
              xNoff_right =  (-0.5 + (1:Noff_right)) * obj.dx_eval./(Nsupersampling);
               offset = [-xNoff_right(end:-1:1),xNoff_right];
           else %insert 0 for odd amount of offset points
              xNoff_right =  (-0.5 + (1:Noff_right)) * obj.dx_eval./(Nsupersampling+1);
               offset = [-xNoff_right(end:-1:1),0,xNoff_right];
           end
           
           temp =  (obj.x_eval*ones(1,Nsupersampling)) +...
                   (ones(Neval,1) * offset);
           
           x_eval_ss  =  reshape (temp, [Nsupersampling*Neval,1]);                       
           %% evaluate data:
           dat_eval_ss = obj.eval (obj_p, data_name, x_eval_ss);
           %% supersampling with Box-Filter:
           temp = reshape(dat_eval_ss, [Neval,Nsupersampling]);           

           dat_eval = mean(temp,2);
        end        
        %%
        function dat_eval = eval(~,obj_p,data_name,x_eval)           
           %cells of the scanning points
           cell_of_xj_eval = cell_structure(obj_p,x_eval);
           %cells of the particles
           cell_of_xj_sph  = cell_structure(obj_p,obj_p.Xj);
           
           %search all neigbhours of xj_eval
           kp=0;   
           pij_eval = []; %[k_eval, k_sph]
           NN = size(x_eval,1);
           for k = 1:NN                
                  cellshift = lookup_cellshift(obj_p,inf);
                  i_neighbours = find(ismember(cell_of_xj_sph,cell_of_xj_eval(k)+cellshift)); 
                  n = length(i_neighbours);
                  pij_eval(kp+(1:n),:) =[k*ones(n,1),i_neighbours];
                  kp=kp+n;
           end
           %compute radius
           xij_h_eval = (x_eval(pij_eval(:,1),:)-obj_p.Xj(pij_eval(:,2),:));
           rij_eval   = sum(xij_h_eval.^2,2).^0.5;
           %compute kernel
           r = rij_eval ./ obj_p.hj(pij_eval(:,2));
           
           I = r < obj_p.kernel_cutoff; 
           
           Wij_eval    = zeros(size(r));
           Wij_eval(I) = obj_p.fw(r(I))./ (obj_p.hj(pij_eval(I,2)).^obj_p.dim);
           
           %work with data
           if strcmp(data_name,'')
               dat = zeros(obj_p.N,1); %e.g. to plot the position
           else
               dat = obj_p.(data_name);
           end
           dat =  dat.*(obj_p.Vj * ones(1,size(dat,2)));
           dat_eval = zeros(NN,size(dat,2));
           %add up everything

           for k = 1:NN                
                II = ismember(pij_eval(:,1),k);
                dat_eval(k,:) = sum((Wij_eval(II)* ones(1,size(dat,2))...
                                    .* dat(pij_eval(II,2),:)),...
                                    1); %sum all rows
           end
        end
        
        %% %%%  Plotting functions %%% %%
        function plot_data(obj,obj_p,A_quantities)
            %% some functions
            function p=plot_scatter(x,dat,mat)
               colo='gbkrmcyw';
               %each material gets his own color
               
               for m = 1:size(mat,1)
                    I = mat(m,1):mat( m,2);
                    p(m)=plot(x(I),dat(I),'o','color',colo(mod(m,length(colo))+1),'MarkerFaceColor','auto'); 
               end                     
            end    
             
            function p=plot_force(data)
                I=data.Icomp;
                F_all = [data.Fj_tot(I), data.Fj_int(I), data.Fj_diss(I)];
                F_leg = ['F_{tot}    ';'F_{int}    ';'F_{diss}   '];
                if any(data.Fj_phy_diss~=0)
                    F_all = [F_all, data.Fj_phy_diss];
                    F_leg = [F_leg; 'F_{phydiss}'];
                end
                if any(data.Fj_ST~=0)
                    F_all = [F_all, data.Fj_ST];
                    F_leg = [F_leg; 'F_{ST}     '];
                end

                p=bar(data.Xj(I),F_all);
                title('forces')
                xlim([data.Omega(1) data.Omega(2)]);
                legend(F_leg);
            end
            
            function p=plot_massflux(data)
                I=data.Icomp;
                p=bar(data.Xj(I),...
                    [data.drhoj(I),data.drhoj_diss(I)]);
                title('massflux')
                xlim([data.Omega(1) data.Omega(2)]);
                legend('drho','drho_{diss}');
            end
            
           
            function scatterPoints=plot_patches(x,y,z,sizeOfCirlce,opacity)
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

            function p=plot_trisurf(x,y,z,max_r)  
                xy   = complex(x,y);
                tri    = delaunay(x,y);
                e    = abs(xy(tri) - xy(circshift(tri,-1,2)));
                goodt = all(e < max_r,2); %dxmedian(e(:)),2);
                t2   = tri(goodt,:);
                p=trisurf(t2,x,y,z);
                %trimesh(t2,x,y,z)
            end
                     
            function [dat_max,p]=plot_field(x,dat)
               p=quiver(x(:,1),...
                         x(:,2),...
                         dat(:,1),...
                         dat(:,2));
               dat_norm = sum(dat.^2,2).^0.5;
               dat_max = max(dat_norm);
            end
           %% start of routine
           if nargin < 3 %use prescribed quantities to plot if not given
               A_quantities = obj.plot_quantity;
           end
            
           if ~isempty(obj.mfigure)
               %do not change the focos (silent update of the figure)
                set(groot,'CurrentFigure',obj.mfigure)
           else
               figure %open new figure
           end
           nplot = length(A_quantities);
           if nplot <=3    
               nxplot = nplot;
               nyplot =1;
           else
               nxplot = ceil(nplot/2);
               nyplot = 2;
           end
           iplot=1;
           x   = obj_p.Xj;
           t   = obj_p.t;
           mat = obj_p.Imaterial_with_boun;
           title_additive = ['; t=',num2str(t,'%10.3e'),'; N= ',num2str(obj_p.N)];
           clf %clear figure
           
           for quantity = A_quantities;
               ab=subplot(nyplot,nxplot,iplot);
               hold on;
               if quantity == 'x'
                   dat_name = '';
                   name = 'position';
                   limaxes = obj.fixaxes.x;
                   style   = obj.plot_style.x;
               elseif quantity == 'p'
                   dat_name = 'pj';
                   name = 'pressure';
                   limaxes = obj.fixaxes.p;
                   style   = obj.plot_style.p;
               elseif quantity == 'd'
                   name = 'density';
                   dat_name = 'rhoj';
                   limaxes = obj.fixaxes.d;                   
                   style   = obj.plot_style.d;
              elseif quantity == 'm'
                   name = 'massflux';
                   dat_name = 'drhoj';
                   limaxes = obj.fixaxes.m;                   
                   style   = obj.plot_style.m;
               elseif quantity == 'v'
                   name = 'velocity';
                   dat_name = 'vj';
                   limaxes = obj.fixaxes.v;
                   style   = obj.plot_style.v;
               elseif quantity == 'f'
                   name = 'F-total'; %todo plot components
                   dat_name = 'Fj_tot';
                   limaxes = obj.fixaxes.f;
                   style   = obj.plot_style.f;
               elseif quantity == 'e'
                   name = 'energy';
                   dat_name = 'ej';
                   limaxes = obj.fixaxes.e;
                   style   = obj.plot_style.e;
               else
                   warning('variable is not implemented for plotting');
                   keyboard
               end

               
               %plot:
               if obj_p.dim == 1   
                   if quantity == 'f'
                       plot_force(obj_p);
                   elseif quantity == 'm'
                       plot_massflux(obj_p);
                   else                            
                       % plot particles
                       p = plot_scatter(x,obj_p.(dat_name),mat);                        
                       pp   =p(1); % array of plot handles                       
                       names={'SPH-data'}; % array of strings for legend
                       
                       % plot smoothed values
                       if ~isempty(obj.x_eval)
                           %evaluate 
                           dat_eval = obj.eval_ss(obj_p,dat_name);

                           p=plot(obj.x_eval,dat_eval,'-');
                           pp=[pp;p];
                           names=[names; 'smoothed values'];
                       end
                       
                       % exact solution:
                       if ~isempty(obj.exact_sol)
                           obj.exact_sol.compute(t,x);
                           p=plot(x,obj.exact_sol.(dat_name),'r.');
                           pp=[pp;p];
                           names=[names; 'analytical solution'];
                       end
                       
                       %initial state of ghost particle   
                       if ~isempty(obj_p.IOdata.([dat_name,'ghost']))
                            plot(x(obj_p.Ighost),obj_p.IOdata.([dat_name,'ghost']),'kx');                       
                       end
                       legend(pp, names);

                       %mark mirror particle
%                        if ~isempty(obj_particles.bc)
%                             plot(obj_particles.Xj(obj_particles.bc.mirrorParticlesj),dat(obj_particles.bc.mirrorParticlesj),'xr');
%                        end
                       
                   end
                   if ~isempty(limaxes)
                     ylim(limaxes);
                   end                   
               elseif obj_p.dim == 2
                   if quantity == 'x'
                      plot_scatter(x(:,1),x(:,2),mat); 
                      axis equal
                   elseif any(quantity == 'pdem') %perssure, density, energy -> scalar
                        if strcmp (style,'trisurf')
                            plot_trisurf(x(:,1),x(:,2),obj_p.(dat_name),2*max(obj_p.hj));
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
                            sizeOfCirlce = obj_p.hj;
                            plot_patches(x(:,1),x(:,2),obj_p.(dat_name),sizeOfCirlce,opacity)                 
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            colorbar 
                            colormap jet
                            axis equal
                        elseif strcmp (style,'plot3')
                            plot3(x(:,1),x(:,2),obj_p.(dat_name),'o');
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                                 zlim(limaxes)
                            end
                            view([0,-1,0]);
                        else
                            error([style, '- plotstyle is not supported']);
                        end
                   elseif any(quantity == 'vf') % velocity, forces -> field
                           dat_max = plot_field(x(obj_p.Iin,:),obj_p.(dat_name)(obj_p.Iin,:));
                           axis equal
                           title_additive= [title_additive,'; max|',quantity,'|=',num2str(dat_max)];
                   else
                       error([obj.plotstyle, '- plotstyle is not supported']);
                   end                   
               end 
               title([name,title_additive]);
               title_additive='';
               %% plot additional info
               plot_geometry(obj,obj_p);
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
                   a = axis;
                   if ~isempty(boun.damping_area)
                       b = boun.damping_area;
                       rectangle('Position',[b(1), a(3), b(2)-b(1), a(4)-a(3)],...
                                 'EdgeColor',[0,1,0.5]);
                   end
                   
                   if strcmp(boun.type,'noflow')
                        plot([boun.p1,boun.p1],a(3:4),'kx:');
                   elseif strcmp (boun.type(1:2),'nr') && obj_p.dim==1
                        plot([obj_p.Xj(boun.kb),obj_p.Xj(boun.kb)],a(3:4),'bx:');
                   else
                       error('No such boundary type to plot');
                   end
                       
                   
                   axis(a);
               end
                
           elseif obj_p.dim == 2
               
               %axis('equal')                   
               %Omega:
               rectangle('Position',[obj_p.Omega(1,1),obj_p.Omega(2,1),...
                   obj_p.Omega(1,2)-obj_p.Omega(1,1),obj_p.Omega(2,2)-obj_p.Omega(2,1)]); %boundary
               hold on
               for boun = obj_p.bc
                   a = axis;
                   plot([boun.p1(1),boun.p2(1)],[boun.p1(2),boun.p2(2)],'kx:');
                   axis(a);
               end
               %% draw connectivity
               if draw_connectivity               
                   colo='gbkrmcyw';
                   AiX=1:500:obj_p.N;
                   for kk = 1:length(AiX)
                       iX= AiX(kk);
                        neigh_iX = [obj_p.pij(obj_p.pij(:,1)==iX,2);obj_p.pij(obj_p.pij(:,2)==iX,1)]; %neighbours of iX
                        for k=1:size(neigh_iX,1)
                           plot([obj_p.Xj(iX,1);obj_p.Xj(neigh_iX(k),1)],...
                                [obj_p.Xj(iX,2);obj_p.Xj(neigh_iX(k),2)],colo(mod(kk,length(colo))+1));
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
        function read_hdf5(~,obj_data,filename,group)
  
            function x=readVariable(x_name,filename,group)
                 x=h5read(filename,[group,'/',x_name]);            
            end
            filename = [filename,'.h5'];
            
            if ~exist(filename,'file')
                error([filename,' does not exist!']);
            end
            
            if isprop(obj_data,'dim')
                dim = obj_data.dim;
            else
                dim = 2; %for lime-sph
            end
            
            %position
            if dim == 1
                obj_data.Xj = readVariable('x',filename,group); 
                obj_data.vj = readVariable('u',filename,group);                
            else
                obj_data.Xj = [readVariable('x',filename,group),...
                               readVariable('y',filename,group)];
                %velocity                       
                obj_data.vj = [readVariable('u',filename,group),...
                               readVariable('v',filename,group)];
            end
                
            %speed of sound
            %c = readVariable('c',filename,time_str);
            obj_data.c0j = readVariable ('c0',filename,group);
            obj_data.cj = readVariable('c',filename,group);

            %density            
            %rho = readVariable('rho',filename,time_str);
            obj_data.rho0j = readVariable('rho0',filename,group);
            obj_data.rhoj = readVariable('rho',filename,group);
     
            %mass
            mj = readVariable('m',filename,group);
            obj_data.Vj = mj ./ obj_data.rhoj;
            
            %energy
            obj_data.ej = readVariable('e',filename,group);
            % some extra information (only available from a simulation
            % output - and not necessary to create a scenario)
            if isprop(obj_data,'hj')
                % smoothing length
                obj_data.hj = readVariable('h',filename,group);
                %pressure
                obj_data.pj = readVariable('p',filename,group);

                %time
                obj_data.t = str2double(group(2:end));                
            end
            
            % take care for the indices
            N = size(obj_data.Xj,1);
            
            if isempty(obj_data.Iin)
                obj_data.Iin = (1:N)';
            end
            
            if isempty(obj_data.Imaterial)
                obj_data.Imaterial = [1,N]; 
            elseif N>max(max(obj_data.Imaterial)) %add material (of boundary)
                obj_data.Imaterial_with_boun = [obj_data.Imaterial;
                                           max(max(obj_data.Imaterial)),N];
            end
           keyboard
          %  Gamma = readVariable('Gamma',filename,time_str);
          %  Gmod = readVariable('Gmod',filename,time_str);
          %  S = readVariable('S',filename,time_str);
          %  Y0 = readVariable('Y0',filename,time_str);
          %  p = readVariable('p',filename,time_str);
          %  phi = readVariable('phi',filename,time_str);
          %  tauXX = readVariable('tauXX',filename,time_str);
          %  tauXY = readVariable('tauXY',filename,time_str);
          %  tauYY = readVariable('tauYY',filename,time_str);

        end
        %%
        function write_hdf5(obj,obj_p)
            filename = [obj.output_name,'.h5'];
            time = obj_p.t;
            function writeVariable(filename,time,name,data)
                group = ['/',num2str(time),'/',name];
                hdf5write(filename, group, data, 'WriteMode', 'append');
            end
            writeVariable(filename,time,'x',obj_p.Xj(:,1));
            writeVariable(filename,time,'u',obj_p.vj(:,1));
            writeVariable(filename,time,'c',obj_p.cj);
            writeVariable(filename,time,'c0',obj_p.c0j);
            writeVariable(filename,time,'rho',obj_p.rhoj);
            writeVariable(filename,time,'rho0',obj_p.rho0j);
            writeVariable(filename,time,'m',obj_p.mj);                        
            writeVariable(filename,time,'h',obj_p.hj);                        
            writeVariable(filename,time,'p',obj_p.pj);                        

            if obj_p.dim == 2
                writeVariable(filename,time,'y',obj_p.Xj(:,2));
                writeVariable(filename,time,'v',obj_p.vj(:,2));                
            end
            
            %% out
            % Gamma, Gmod, O, S, Vol, Y0, c, c0, e, eh, epsdotXX, epsdotXY,
            % et, h , m, nn, p ,phi, rho, rho0, rhoh, rhot, rotdotXY,
            % stationary, t, tauXX, tauXXpl, tauXXrho2, tauXY, tauXYpl,
            % tauXYrho2, tauYY, tauYYpl,tauYYrho2, taudotXX, taudotXY, 
            % taudotYY, tmp, u, uh, ups, ut, v, vh, vt, x ,y , 
        end
    end

    
end

