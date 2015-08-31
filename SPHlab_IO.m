% SPH for HVI  - Markus Ganser - TU/e - 2015
classdef SPHlab_IO < handle
    % IO classes for the sph-simulation tool
    % most configuration are done in SPHlab_scenario or manually, see
    % constructor for the variables/adjustments handed over 
    
    % - initialize | do | finalize
    %   main functions for plotting and record of conservation properties
    % - save | plot_conservation
    %   record conservation of mass/momentum/energy, plot only if changes
    %   are present
    % - plot | save_pointdata
    %   record the state variables of one particle - to be adjusted in the
    %   function itself
    % - dat_eval
    %   evaluates on an arbitrary point with the SPH discretization schemes
    %   supersampling possible (but not necessary)
    % - plot_trisurf
    %   seperated function to be applicable from the outside
    % - plot_data
    %   main plotting routine with different configurations (see
    %   SPHlab_scenario constructor for full information)
    %   with plot_data("obj_p","var"), a plot with the variable "var" (e.g. p
    %   for pressure) based on the data "obj_p" (the particle class) can be
    %   plottet
    % - draw_geometry:
    %   adds additional information to the plot (e.g. the domain,
    %   cellstructure,..), see flags in this function
    % - draw_points
    %   mark points with a cross, Ind denotes the Indexarray of the
    %   particles to be plottet
    % - savefigure
    %   save a figure in a eps/png/... file
    % - read | write_hdf5
    %   read and writes hdf5 files, compatible with LimeSPH
    
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
        plotconfig %struct
        %
        
        Sfont
        Sint
        
        save_as_figure
        figure_format
        %movie
        save_as_movie
        movie_name
        vidObj
       
        %conservation variables
        con_dt
        con_mass
        con_momentum
        con_ang_momentum
        con_momentum_diss
        con_e_kin
        con_e_kin_diss
        con_e_pot
        con_e_pot_diss
        
        %some internal flags/variables
        t_last_plot
        t_last_save
     
        %data of (one) point (struct)
        pointdata
        
        %exact solution
        exact_sol
    end
    
    methods
        
        %% %%% constructor
        function obj = SPHlab_IO(obj_scen)
            if nargin == 1 %else, use only the functions
              %% some IO properties
              obj.save_as_movie = obj_scen.save_as_movie;
              obj.save_as_figure = obj_scen.save_as_figure;
              obj.figure_format = obj_scen.figure_format;
              obj.movie_name    = get_movie_name(obj_scen);
              obj.save_dt       = obj_scen.save_dt;
              obj.plot_dt       = obj_scen.plot_dt;
              obj.plot_style    = obj_scen.plot_style;
              obj.plot_quantity = obj_scen.plot_quantity;                                          
              obj.exact_sol     = obj_scen.exact_sol;                                          
              obj.plotconfig    = obj_scen.plotconfig;
              obj.fixaxes       = obj_scen.fixaxes;
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
            end 
          %% conservation variables              
          obj.con_mass = [];
          obj.con_momentum =[];
          obj.con_ang_momentum =[];
          obj.con_momentum_diss =[];
          obj.con_dt = [];
          obj.con_e_kin =[];
          obj.con_e_kin_diss =[];
          obj.con_e_pot =[];
          obj.con_e_pot_diss =[];

          %%
          obj.t_last_plot = -inf;
          obj.t_last_save = -inf;
          obj.pointdata=struct('t',[],'hj',[],'Omegaj',[],'rhoj',[],'pj',[],'vj',[],'ej',[],'drhoj_tot',[],'drhoj_diss',[],'maxp',[]);  

          %for latex plots:
            %----------------
           FontSize=20;%15;
           obj.Sfont= struct('FontUnits','points',...
            'FontSize',FontSize,...
            'FontName','Times');
           obj.Sint = struct('interpreter','latex');
            %------------------
        end
        
        %% %%% general functions
        function initialize(obj)

            %% plot
            if ~strcmp(obj.plot_quantity,'')  %plot only, if plotstlye is specified
                 obj.mfigure = figure('Color',[1 1 1]);
                 set(gca,'DataAspectRatio',[1,1,1]);
                 if strcmp(get(0,'DefaultFigureWindowStyle'),'normal') %only when not docked
                     if ~isempty(obj.plotconfig.figuresize)
                         obj.mfigure.Units = 'centimeters';
                         obj.mfigure.Position =obj.plotconfig.figuresize;                         
                         obj.mfigure.PaperPosition =obj.plotconfig.figuresize;
                     else
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
                            figpos(2)=0.1;
                            figpos(4)=0.5;
                            obj.mfigure.Position =figpos;
                        end
                     end
                 end
            end

            
            %% movie
            if obj.save_as_movie
%                 set(gca,'DataAspectRatio',[1,1,1]);
%                 obj.mfigure.Units = 'normalized';
%                 obj.mfigure.Position = [0 0 1 1];
                obj.vidObj  = VideoWriter(obj.movie_name,'MPEG-4');
                obj.vidObj.FrameRate = 5;

                open(obj.vidObj);                
            end   
            
            %% output
           

        end
        %%
        function do (obj, obj_particles,force_do,flag_conservation)
            % plot and save data from obj_particles
            % force_do: plot and save data independent of the configuration 
            % flag_conservation: save conservation of mass/momentum/energy                        
            
            if nargin <3
                force_do = false;
            end
            if nargin < 4
                flag_conservation = true;
            end
            %% plotting
            if (force_do||(((obj_particles.t-obj.t_last_plot) > obj.plot_dt) ...
                    && ~isempty(obj.plot_quantity)))
                plot_data (obj,obj_particles);
                if obj.save_as_movie
                    currFrame = getframe(obj.mfigure);
                    writeVideo(obj.vidObj,currFrame);
                end   
                if obj.save_as_figure
                    obj.savefigure(obj_particles.t)
                end
                obj.t_last_plot = obj_particles.t;
            end  
            
            %% write output
            if (force_do||((obj_particles.t-obj.t_last_save) > obj.save_dt)) ...
                            && obj.write_data
                 write_hdf5(obj,obj_particles)
                    %name = ['data/out',num2str(obj_particles.dt),'.mat'];               
                    % save(name, 'obj_particles');   
                 obj.t_last_save = obj_particles.t;
            end
            
            %check for conservation of mass and momentum            
            if flag_conservation
                save_conservation (obj,obj_particles);
            end
            % save data from one special point
            save_pointdata(obj,obj_particles)

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
%             disp ('some IO plots turned off...');
           plot_conservation(obj)
           plot_pointdata(obj);
            

            
        end
        %%
        function save_conservation (obj,obj_p)
            % config:
            A=obj_p.Iin; %which particles shall be considered
            
            %time
            obj.con_dt = [obj.con_dt;
                          obj_p.t];
            
            % conservation of mass
            massj = obj_p.rhoj(obj_p.Iin).*obj_p.Vj(obj_p.Iin);
            obj.con_mass = [obj.con_mass;...
                            sum(massj(A))];
                        
            %% conservation of momentum (particle inside)
            momj  = massj * ones(1,obj_p.dim) .* (obj_p.vj_half(obj_p.Iin,:));
            obj.con_momentum = [obj.con_momentum;...
                                sum(momj(A,:),1)];
                            
            %% absorbed momentum:                          
            nrm  = false;
            Inrc = 0;
            for boun=obj_p.bc
                if strcmp(boun.type,'nrm')
                   nrm = true;
                elseif  strcmp(boun.type,'nrc')                   
                   % nrc:     
                   dI = boun.dvj_diss.*(obj_p.mj(boun.kb)*ones(1,obj_p.dim));
                   Inrc = Inrc + sum(dI,1) *obj_p.dt;   
                end
            end
            
            if (nrm && obj_p.compGhost) 
                % nrm: dissipation of momentum  (at half step)          
                momj_diss_add   = obj_p.IOdata_h.vjghost.*(obj_p.mj(obj_p.Ighost)*ones(1,obj_p.dim));                
                momj_diss_after = obj_p.vj_half(obj_p.Ighost,:).*(obj_p.mj(obj_p.Ighost)*ones(1,obj_p.dim));
                Inrm = sum(momj_diss_after-momj_diss_add,1);
            else
                Inrm = 0;
            end
            
            if isempty(obj.con_momentum_diss) %time integration
                previous = zeros(1,obj_p.dim); 
            else
                previous = obj.con_momentum_diss(end,:);
            end
            obj.con_momentum_diss =[obj.con_momentum_diss;...
                                    Inrc+Inrm+previous];                      
            
            %% angular momentum
            if obj_p.dim == 2
                %center of mass:
%                 Xj_mass_center = [sum(obj_p.mj.*obj_p.Xj(obj_p.Iin,1))/sum(obj_p.mj),...
%                     sum(obj_p.mj.*obj_p.Xj(obj_p.Iin,2))/sum(obj_p.mj)]                
                % angular momentum:
                amomj = obj_p.mj(obj_p.Iin).*((obj_p.Xj(obj_p.Iin,1)).*obj_p.vj_half(obj_p.Iin,2)...
                       -(obj_p.Xj(obj_p.Iin,2)).*obj_p.vj_half(obj_p.Iin,1));
                obj.con_ang_momentum = [obj.con_ang_momentum;...
                                         sum(amomj)];
            end
            
            %energy:
            if obj_p.dim==2
                norm_vj_sq = (obj_p.vj_half(obj_p.Iin,1).^2 + obj_p.vj_half(obj_p.Iin,2).^2);
            else
                norm_vj_sq = obj_p.vj_half(obj_p.Iin,1).^2;
            end
            e_kin = 0.5 .* obj_p.mj(obj_p.Iin).* norm_vj_sq;
            e_pot = obj_p.mj(obj_p.Iin).*obj_p.ej_half(obj_p.Iin);
            obj.con_e_kin =[obj.con_e_kin;
                            sum(e_kin(A))];
            obj.con_e_pot =[obj.con_e_pot;
                            sum(e_pot(A))];
                        
            % energy dissipation on the boundary (absorbed energy)
            
            %nrm:
            if (nrm && obj_p.compGhost)
                %dv = data.vj(data.Ighost) - data.IOdata.vjghost;
                if obj_p.dim==2
                    norm_vj_sq      = obj_p.vj_half(obj_p.Ighost,1).^2 +...
                                      obj_p.vj_half(obj_p.Ighost,2).^2;
                    norm_vj_sq_prev = obj_p.IOdata_h.vjghost(:,1).^2 + ...
                                      obj_p.IOdata_h.vjghost(:,2).^2;
                else
                    norm_vj_sq      = obj_p.vj_half(obj_p.Ighost,1).^2;
                    norm_vj_sq_prev = obj_p.IOdata_h.vjghost(:,1).^2;
                end
                
                e_kin_ghost      = 0.5 .* obj_p.mj(obj_p.Ighost).* norm_vj_sq;
                e_kin_ghost_prev = 0.5 .* obj_p.mj(obj_p.Ighost).* norm_vj_sq_prev;                
                e_kin_ghost_diff = e_kin_ghost - e_kin_ghost_prev;
                
                e_pot_ghost      = obj_p.ej(obj_p.Ighost).*obj_p.mj(obj_p.Ighost);
                e_pot_ghost_prev = obj_p.IOdata_h.ejghost.*obj_p.mj(obj_p.Ighost);
                e_pot_ghost_diff = e_pot_ghost - e_pot_ghost_prev;                
               
                Ekin_nrm = sum(e_kin_ghost_diff);
                Epot_nrm = sum(e_pot_ghost_diff);
            else
                Ekin_nrm = 0;
                Epot_nrm = 0;
            end
            
            
            %nrc:
            Ekin_nrc=0;
            Epot_nrc=0;
            for boun=obj_p.bc
                % nrc:     
                if  strcmp(boun.type,'nrc')                    
                   %kinetic energy:
                   dEkin    = obj_p.vj(boun.kb,:).*boun.dvj_diss.*(obj_p.mj(boun.kb)*ones(1,obj_p.dim));
                   Ekin_nrc = Ekin_nrc + sum(sum(dEkin.^2,2).^0.5) *obj_p.dt;  
                   %potential energy:
                   dEpot    = boun.dej_diss .*obj_p.mj(boun.kb);
                   Epot_nrc = Epot_nrc + sum(dEpot) *obj_p.dt;
                   
                end
            end
            
            % save data
            if isempty(obj.con_e_kin_diss)
                    previous_kin = 0;
                    previous_pot = 0;
                else                
                    previous_kin = obj.con_e_kin_diss(end);
                    previous_pot = obj.con_e_pot_diss(end);
            end
            
            obj.con_e_kin_diss =[obj.con_e_kin_diss;
                                     Ekin_nrc + Ekin_nrm + previous_kin];
            obj.con_e_pot_diss =[obj.con_e_pot_diss;
                                     Epot_nrc + Epot_nrm + previous_pot];

        end
        %%
        function plot_conservation (obj)
            if isempty(obj.con_dt)
                disp('no conservation record available')
                return
            end
            
            %set flags
            plotmass = any(diff(obj.con_mass) ~= 0);
            plotmomentum = ~isempty(obj.con_momentum);
            plotangmomentum = ~isempty(obj.con_ang_momentum);
            plotenergy = ~isempty(obj.con_e_pot) && ~isempty(obj.con_e_kin);
            
            %subplot counter:
            fig=figure;
            nplot = plotmass + plotmomentum + plotangmomentum + plotenergy;            
            iplot = 1;
            
            %% -------------------------------------
            if plotmass  % plot only if change of mass occurs              
                subplot(nplot,1,iplot)
                box on;
                iplot=iplot+1;
                plot(obj.con_dt,obj.con_mass)
                h1=xlabel('$t$');
                h2=ylabel('$m$');
                if obj.plotconfig.latexplot
                   set(h1,obj.Sfont,obj.Sint);
                   set(h2,obj.Sfont,obj.Sint);
                   % axis
                   set(gca,obj.Sfont);
                else
                   title('evolution of mass');
                end
            end
            %% -------------------------------------
            %momentum
            if plotmomentum
                subplot (nplot,1,iplot)
                grid on;
                hold on;
                box on;
                iplot=iplot+1;
                %components (in 2d)
                dim = size(obj.con_momentum,2);            
                if dim ==1
                    plot(obj.con_dt,obj.con_momentum); 
                    leg = '$I_{int}$   ';
                    %dissipation:
                    if ~isempty(obj.con_momentum_diss)
                        plot(obj.con_dt,obj.con_momentum_diss,':');   
                        plot(obj.con_dt,obj.con_momentum_diss+obj.con_momentum,'k-.');   
                        leg = [leg;'$I_{abs}$   ';'$\sum I$    '];
                    end
                    mom_norm = abs(obj.con_momentum);           

                else
                    mom_norm = sum(obj.con_momentum.*obj.con_momentum,2).^0.5;
                    plot(obj.con_dt,mom_norm);    
                    leg = '$|I|$       ';
                    %components:
                    for i = 1: dim
                        plot(obj.con_dt,obj.con_momentum(:,i)) %x
                        leg = [leg;'I_{x',num2str(i),'}      '];
                    end
                    %dissipation:
                    mom_norm_diss = sum(obj.con_momentum_diss.*obj.con_momentum_diss,2).^0.5;
                    if ~isempty(mom_norm_diss)
                        plot(obj.con_dt,mom_norm_diss,':');   
                        plot(obj.con_dt,mom_norm_diss+mom_norm,'k-.');   
                        leg = [leg;'$|I_{abs}|$ ';'$|\sum I|$  '];
                    end                                        
                end
                h1=xlabel('$t$'); 
                h2=ylabel('$I$');
                h3=legend(leg);
                if obj.plotconfig.latexplot
                   set(h1,obj.Sfont,obj.Sint);
                   set(h2,obj.Sfont,obj.Sint);
                   set(h3,obj.Sfont,obj.Sint);
                   % axis
                   set(gca,obj.Sfont);
                else
                   title('evolution of momentum (at half step)');
                end

            end
            %% -------------------------------------
            if plotangmomentum
                subplot(nplot,1,iplot)
                box on;
                grid on;
                hold on;
                iplot=iplot+1; 
                plot(obj.con_dt,obj.con_ang_momentum)
                h1=xlabel('t'); 
                h2=ylabel('angular momentum');
                if obj.plotconfig.latexplot
                   set(h1,obj.Sfont,obj.Sint);
                   set(h2,obj.Sfont,obj.Sint);
                   set(h3,obj.Sfont,obj.Sint);
                   % axis
                   set(gca,obj.Sfont);
                else
                   title('evolution of angular momentum');
                end
                ang_mom_norm = sum(obj.con_ang_momentum.*obj.con_ang_momentum,2).^0.5;
            else
                ang_mom_norm=[];
            end
            %% -------------------------------------
            if plotenergy
                %energy
                subplot(nplot,1,iplot)   
                hold on;
                grid on
                box on;
                sum_e  = obj.con_e_pot+obj.con_e_kin;
                etot  = obj.con_e_pot+obj.con_e_kin;
                plot(obj.con_dt,obj.con_e_pot,...
                     obj.con_dt,obj.con_e_kin);                     
                 
                 leg = ['$e_{pot}$   ';'$e_{kin}$   '];
                 if ~isempty(obj.con_e_pot_diss)
                     plot(obj.con_dt,obj.con_e_pot_diss,'--',...
                        obj.con_dt,obj.con_e_kin_diss,'--');          
                     sum_e= sum_e+(obj.con_e_pot_diss+obj.con_e_kin_diss); 
                     leg = [leg;'$e_{pot,ab}$';'$e_{kin,ab}$'];
                 end
                 plot(obj.con_dt,sum_e,'k-.');
                 leg = [leg;'$\sum e$    '];

                 h1=xlabel('$t$'); 
                 h2=ylabel('$e$');
                 h3=legend(leg);
                 if obj.plotconfig.latexplot
                   set(h1,obj.Sfont,obj.Sint);
                   set(h2,obj.Sfont,obj.Sint);
                   set(h3,obj.Sfont,obj.Sint);
                   % axis
                   set(gca,obj.Sfont);
                 else
                  title('energy (at half step)');
                 end
            end
            %% -------------------------------------
            %move figure to the left side (only when not docked)
            if strcmp(get(0,'DefaultFigureWindowStyle'),'normal')
                figpos=fig.Position;
                figpos(2)=0;
                fig.Position =figpos;
            end
          
            % output in console:
            if ~isempty(mom_norm)
                disp (['## rel. diss. of momentum = ',...
                    num2str(abs((mom_norm(end)-mom_norm(1))/mom_norm(1))),' ##']);                
            end
            if ~isempty(ang_mom_norm)
                disp (['## rel. diss. of angular momentum = ',...
                        num2str(abs((ang_mom_norm(end)-ang_mom_norm(1))/ang_mom_norm(1))),' ##']);
            end
            if ~isempty(etot)
                disp (['## rel. diss. of tot. energy = ',...
                        num2str(abs((etot(end)-etot(1))/etot(1))),' ##']);
            end
        end
        %%
        function save_pointdata(obj,obj_particles)
             %save data of particle iparticle
            iparticle=289;%100;%272;%;251;
            
            if iparticle > obj_particles.N
                return
            end
            obj.pointdata.t(end+1)= obj_particles.t;
            if length(obj_particles.hj)>1
                 obj.pointdata.hj(end+1)= obj_particles.hj(iparticle);                
            else
                 obj.pointdata.hj(end+1)= obj_particles.hj(1);                
            end
            if length(obj_particles.Omegaj)>1
                 obj.pointdata.Omegaj(end+1)= obj_particles.Omegaj(iparticle);                
            else
                 obj.pointdata.Omegaj(end+1)= obj_particles.Omegaj(1);                
            end
            obj.pointdata.rhoj(end+1)= obj_particles.rhoj(iparticle);
            obj.pointdata.pj(end+1)= obj_particles.pj(iparticle);
            obj.pointdata.vj(end+1)= obj_particles.vj(iparticle);
            obj.pointdata.drhoj_diss(end+1)= obj_particles.drhoj_diss(iparticle);
            obj.pointdata.drhoj_tot(end+1)= obj_particles.drhoj_tot(iparticle);
            obj.pointdata.ej(end+1)= obj_particles.ej(iparticle);           
            obj.pointdata.maxp(end+1)=max(abs(obj_particles.pj));
        end
        %%         
        function plot_pointdata(obj)
           p = struct('pq',[]);
            %--------------------
           %config:
           p(1).pq='pj';%'rhoj';
%            p(2).pq='vj';
        
           marker='-';
           %----------------------                
           if isempty(obj.pointdata.t)                                
                return               
           end
           np=size(p,2);
           figure(30)
             % plot pointdata
            
            for ip=1:np
                subplot(1,np,ip)
                hold on;
                plot(obj.pointdata.t,obj.pointdata.(p(ip).pq),marker);            hold on;
                xlabel('t'),ylabel(p(ip).pq);            
            end
%            figure(31)
%            plot(obj.pointdata.t,obj.pointdata.maxp,'-b');
%            hold on;
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
           if isempty(pij_eval)
              error('no particles in the region of evalution'); 
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

           for k = 1:NN   %evaluate point for point             %Violeau (5.86)
                II = ismember(pij_eval(:,1),k);
                if sum(II) == 0
                    continue
                end
                %normalization constant
                sumW = sum(obj_p.Vj(pij_eval(II,2)).*Wij_eval(II));
                if sumW == 0 % no data to sum up, so just let it be zero
                    continue              
                end
                dat_eval(k,:) = sum((Wij_eval(II)* ones(1,size(dat,2))...
                                    .* dat(pij_eval(II,2),:)),...
                                    1)/sumW; %sum all rows
           end
        end        
        %% %%%  Plotting functions %%% %%        
        function p = plot_trisurf(~,x,y,z,max_r)  
            %seperate because it is needed as an extra routine for a
            %comparison analysis, too.
            xy    = complex(x,y);
            tri   = delaunay(x,y);
            e     = abs(xy(tri) - xy(circshift(tri,-1,2)));
            goodt = all(e < max_r,2); %dxmedian(e(:)),2);
            t2    = tri(goodt,:);
            p     = trisurf(t2,x,y,z,'EdgeColor','none');            
            %trimesh(t2,x,y,z)
        end  

        %% main plotting function:
        function plot_data(obj,obj_p,A_quantities,flag_drawnow,flag_cla)
            %% some functions
            function p=plot_scatter(x,dat,mat,marker)
%                colo='gkbkbrmcy';
               colo='kbrbrb';
               %each material gets his own color
               
               for m = 1:size(mat,1)
                    I = mat(m,1):mat( m,2);
                    p(m)=plot(x(I),dat(I),marker,'color',colo(mod(m,length(colo))+1),'MarkerFaceColor','auto'); 
               end                     
            end    
             
            function p = plot_force(data)
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

                p = bar(data.Xj(I),F_all);
                title('forces');
                xlim([data.Omega(1) data.Omega(2)]);
                legend(F_leg);
            end
            
            function p = plot_massflux(data)
                I=data.Icomp;
                p=bar(data.Xj(I),...
                    [data.drhoj_tot(I),data.drhoj_diss(I)]);
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
                     
            function [dat_max,p]=plot_field(x,dat)
               p=quiver(x(:,1),...
                         x(:,2),...
                         dat(:,1),...
                         dat(:,2));
               dat_norm = sum(dat.^2,2).^0.5;
               dat_max = max(dat_norm);
            end
            
            function p=plot_torus(Ar,Az,Ad)
                warning('super slow')
                keyboard
                N = size(Ad,1);
                for i=1:10:N
                     d = Ad(i);
                     r  = Ar(i);
                     zz = Az(i);
                     [u,v]=meshgrid(0:30:60);                 
                     X=(r+d*cosd(v)).*cosd(u);
                     Y=(r+d*cosd(v)).*sind(u);
                     Z=d*sind(v)+zz;
                     p=surfc(X,Y,Z);
                     hold on
                     disp([num2str(i),'/',num2str(N)]);
                end
            end
            
            function p=plot_ring(xx,yy)
                N = size(xx,1);
                Ntheta=20;
                theta = linspace(0,pi/3,Ntheta);
                for i=1:N 
                    r = abs(yy(i));
                    Y = r*(cos(theta));
                    Z = r*(sin(theta));
                    X = xx(i)*ones(Ntheta,1);
                    p=plot3(X,Y,Z,'-b');     
                    disp([num2str(i),'/',num2str(N)]);
                end
            end %2d-rings in the 3dspace for axisymmetric plots
            
            function p=plot_ringcloud(xx,yy,mat)
               colo='gbkrmcy';
               %each material gets his own color
               Ntheta=20; %amount of points in the axis
               Theta_minus = 0;
               Theta_plus  = pi/2;
               for m = 1:size(mat,1)
                    I = mat(m,1):mat( m,2);                
                    N=size(xx(I),1);
                    theta = linspace(Theta_minus,Theta_plus,Ntheta);
                    r = abs(yy(I));
                    Y = reshape(r*(cos(theta)),1,N*Ntheta);
                    Z = reshape(r*(sin(theta)),1,N*Ntheta);
                    X = reshape(xx(I)*ones(1,Ntheta),1,N*Ntheta);
                    p = plot3(X,Y,Z,'.','color',colo(mod(m,length(colo))+1));%,'MarkerFaceColor','auto');  
               end
            end            
            
             function add_fluxes(obj_p,datname)               
               %for boun
               datnamedisp=['d',datname,'_{diss}'];
               scalefactor = 10;
               dx_factor = 1;
               len=(obj_p.Omega(1,2)-obj_p.Omega(1,1))/scalefactor;
               for iboun = 1:size(obj_p.bc,2)                   
                   boun=obj_p.bc(iboun);
                   if ~strcmp(boun.type,'nrc') %if not a nr-bc continue
                       continue
                   end
                   if ~isfield(boun,'kb') %handling older data sets (just do not plot that)
                       continue
                   end
                   bp     = boun.bp;
                   normal = boun.outer_normal;
                   kbb    = find(boun.kb);
                   X      = obj_p.Xj(kbb,:);
                   NN     = size(X,1);               
                   dx     = dx_factor *mean(obj_p.Vj(kbb).^(1/obj_p.dim));
                   if obj_p.dim==1
                       ax = axis;
                       y0 = 0.5*(ax(4)+ax(3));
                       dat = sum(boun.(['d',datname,'_diss']));%/5e8;                                            
                         
                     %plotting:                  
                       xright = mean(X)+(ax(2)-ax(1))/scalefactor*normal;
                       plot(xright*[1,1],[y0,y0+dat],'r-x')
                       text(xright*1.05,y0,[datnamedisp,'=',num2str(dat,'%10.2e')])
                       axis(ax)
                   else %2d                             
                      if normal(1)==0  % discretisize vectical 
                           kx=1;
                           ky=2;                       
                       elseif normal(2)==0
                           kx=2;
                           ky=1;
                       else
                           error('todo');
                       end
                       xleft  = min(X(:,kx))-dx/2;
                       xright = max(X(:,kx))+dx/2;
                       Nc = ceil((xright-xleft)/dx);
    %                    dx =(xright-xleft)/n;

                       cell_of_X = floor(ones(NN,1)*(Nc ./ (xright-xleft))'...
                        .* (X(:,kx) - ones(NN,1)*(xleft)'))...
                        +1; %in which sector is the particle

                       %compute data (simply sum up all fluxes in the slide
                       dat=zeros(Nc,obj_p.dim);
                       for i = 1:Nc
                          dat(i,:) = sum(boun.(['d',datname,'_diss'])(cell_of_X==i,:),1);%/5e8;                                            
                       end

                       %just take the norm
                       dat = sum(dat.^2,2).^0.5;
                       %compute max and scale data
                       datmax=(max(abs(dat)));
                       if datmax == 0
                           continue
                       end
                       y0= bp(ky)+normal(ky)*1.5*len;                         
                       %plot reference box
                       if kx==1 %to the top/bottom
                           rectangle('Position',[xleft,y0-len,xright-xleft,2*len],'LineStyle',':')
                           plot([xleft;xright],[y0,y0],'b');
                       else %to the right/left
                           rectangle('Position',[y0-len,xleft,2*len,xright-xleft],'LineStyle',':')
                           plot([y0,y0],[xleft;xright],'b');
                       end
                       %normalize
                       dat1= normal(ky)*dat*len/datmax;
                       %plot only non-zero values (espacially when no
                       %particle belongs to an Nc-cell
                       xxx=xleft:dx:xright;
                       I=dat1~=0;
                       %plot data
                       if kx==1 %to the top/bottom
                            plot(xxx(I),y0+dat1(I)','rx')
                            %text(xright,y0,['max ',datnamedisp,'=',num2str(datmax,'%10.2e')])
                       else %to the right/left
                            plot(y0+dat1(I)',xxx(I),'rx')
                            %text(y0,xright+0.2*len,['max ',datnamedisp,'=',num2str(datmax,'%10.2e')])

                       end 
                   end
               end   
             end
             
          %% --------------------------------------------
          %% --------------------------------------------  
          %% start of routine
          if nargin < 3 %use prescribed quantities to plot if not given
               A_quantities = obj.plot_quantity;
          end
           %default: dont wait for drawing the figure and clear axes
           if nargin < 4
               flag_drawnow = true;
           end           
           if nargin < 5
               flag_cla = true;
           end
            
           if ~isempty(obj.mfigure)
               %do not change the focus (silent update of the figure)
              %  set(groot,'CurrentFigure',obj.mfigure)               
           else
               figure('Color',[1 1 1]); %open new figure
           end
           
           nplot = length(A_quantities);
           if nplot <=3    
               nxplot = nplot;
               nyplot = 1;
           else
               nxplot = ceil(nplot/2);
               nyplot = 2;               
           end
           
           %transpose plots
           if obj.plotconfig.transpose
               temp   = nxplot;
               nxplot = nyplot;
               nyplot = temp;
           end
           %---------------
           
           iplot=1;
           x   = obj_p.Xj;
           t   = obj_p.t;
           mat = obj_p.Imaterial_with_boun;
           title_additive = ['; t=',num2str(t,'%10.3e'),'; N= ',num2str(obj_p.N)];
           
           
           for quantity = A_quantities;
               if nplot > 1
                  subplot(nyplot,nxplot,iplot);
               end
               if flag_cla
                  cla %clear axis
               end
               
               hold on;
               axis normal
               if quantity == 'x'
                   dat_name = '';
                   name = 'position';
                   axisname = '';
                   limaxes = obj.fixaxes.x;
                   style   = obj.plot_style.x;
               elseif quantity == 'p'
                   dat_name = 'pj';
                   name = 'pressure';
                   axisname = 'p';%'pressure $[Pa]$';
                   limaxes = obj.fixaxes.p;
                   style   = obj.plot_style.p;
               elseif quantity == 'd'
                   name = 'density';
                   if strcmp(obj_p.scheme,'a') 
                        dat_name = 'rhoj_real'; %plot real density in the axis-symetric case
                   else
                       dat_name = 'rhoj';
                   end
                   axisname = '$\rho$';
                   limaxes = obj.fixaxes.d;                   
                   style   = obj.plot_style.d;
              elseif quantity == 'm'
                   name = 'massflux';
                   dat_name = 'drhoj_tot';
                   axisname = '$\dot{\rho}$';
                   limaxes = obj.fixaxes.m;                   
                   style   = obj.plot_style.m;
               elseif (quantity == 'v') || (quantity == 'u')
                   name = 'velocity';
                   dat_name = 'vj';
                   axisname = 'u';%'velocity $[m/s]$';
                   limaxes = obj.fixaxes.v;
                   style   = obj.plot_style.v;
               elseif quantity == 'f'
                   name = 'F-total'; %todo plot components
                   dat_name = 'Fj_tot';
                   axisname = '$F$';
                   limaxes = obj.fixaxes.f;
                   style   = obj.plot_style.f;
               elseif quantity == 'e'
                   name = 'energy';
                   dat_name = 'ej';
                   axisname = '$e$';
                   limaxes = obj.fixaxes.e;
                   style   = obj.plot_style.e;
               elseif quantity == 'c'
                   name = 'speed of sound';
                   dat_name = 'cj';
                   axisname = '$c$';
                   limaxes = obj.fixaxes.c;
                   style   = obj.plot_style.c;
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
                       marker = '.';
                       p = plot_scatter(x,obj_p.(dat_name),mat,marker);                        
                       pp = p(1); % array of plot handles                       
                       names={'SPH-particles'}; % array of strings for legend
                       % exact solution:
                       if ~isempty(obj.exact_sol)
                           xx=sort(x);
                           obj.exact_sol.compute(t,xx);
                           p=plot(xx,obj.exact_sol.(dat_name),'r-');                           
                           pp=[pp;p];
                           names=[names; 'analytical solution'];
                       end
                       
                       % plot smoothed values
                       if ~isempty(obj.x_eval)
                           %evaluate 
                           dat_eval = obj.eval_ss(obj_p,dat_name);

                           p=plot(obj.x_eval,dat_eval,'-');
                           pp=[pp;p];
                           names=[names; 'smoothed values'];
                       end
                       

                       
                       %initial state of ghost particle   
                       if ~isempty(obj_p.IOdata.([dat_name,'ghost']))
                            plot(x(obj_p.Ighost),obj_p.IOdata.([dat_name,'ghost']),'kx');                       
                       end
                       
                       %plot legend only if more than one property is in
                       %the plot
                       if size(names,1)>1
                           if quantity=='d'
                               h=legend(pp, names);
                               set(h,obj.Sfont,obj.Sint);
                           end
                       end

                       %mark mirror particle
%                        if ~isempty(obj_particles.bc)
%                             plot(obj_particles.Xj(obj_particles.bc.mirrorParticlesj),dat(obj_particles.bc.mirrorParticlesj),'xr');
%                        end
                   end
                   if ~isempty(limaxes)
                     ylim(limaxes);
                   else
                       ylim([-inf,inf]);
                   end  
                   
                   h1=xlabel('$x$');
                   h2=ylabel(axisname);
                   if obj.plotconfig.latexplot
                     set(h1,obj.Sfont,obj.Sint);
                     set(h2,obj.Sfont,obj.Sint);
                       % axis
                     set(gca,obj.Sfont);
                   end
                   
               elseif obj_p.dim == 2
                   addcolorbar = false;
                   if quantity == 'x'
                      if strcmp (style,'ring')         
                            plot_ring(x(:,1),x(:,2));
                            view([-0.3,-1,1]);                                
                      elseif strcmp (style,'ringcloud')         
                            plot_ringcloud(x(:,1),x(:,2),mat);
                            view([-0.3,-1,1]);   
%                             daspect([1,1,1])
                      elseif strcmp (style,'torus')         
                            sizeOfCirlce = obj_p.hj;
                            plot_torus(x(:,1),x(:,2),sizeOfCirlce);
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            view([-0.3,-1,1]);                                
                      else %default
                          if flag_cla %default marker                              
                              marker = '.';
                          else        %alternative marker (in order to compare two solutions)
                              marker ='x';
                          end
                          plot_scatter(x(:,1),x(:,2),mat,marker); 
                      end
                      axis equal;
                   elseif size(obj_p.(dat_name),2)==1 %pressure, density, energy, massflux, speed of sound -> scalar
                        if strcmp (style(1:5),'trisu')
                            if strcmp(dat_name,'rhoj') && strcmp(obj_p.scheme,'a')
                                %scale density
                                plot_trisurf(obj,x(:,1),x(:,2),...
                                    obj_p.(dat_name)./(2*pi*abs(x(:,2))),...
                                    2*max(obj_p.hj));
                            else %standard
                                plot_trisurf(obj,x(:,1),x(:,2),obj_p.(dat_name),2*max(obj_p.hj));
                            end
                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            if length(style)>7 %make it from a perspective, e.g. 'trisurf3'
                                view([1,0.5,0.3]);
                            else     
                                view(2)
                                daspect([1,1,1])
                            end
                            
                            shading interp
%                             if ~isempty(limaxes)
%                                  zlim(limaxes)
%                             end
                            addcolorbar = true;
                    
                        elseif strcmp (style,'patches')                            	
                            opacity = 0.3;
                            sizeOfCirlce = obj_p.hj;

                            plot_patches(x(:,1),x(:,2),obj_p.(dat_name),sizeOfCirlce,opacity);

                            if ~isempty(limaxes)
                                 caxis(limaxes)
                            end
                            addcolorbar = true;
                            daspect([1,1,1])
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
                   elseif size(obj_p.(dat_name),2)==2  % velocity, forces -> field
                           dat_max = plot_field(x(obj_p.Iin,:),obj_p.(dat_name)(obj_p.Iin,:));
                           daspect([1,1,1])
                           title_additive= [title_additive,'; max|',quantity,'|=',num2str(dat_max)];
                   else
                       error([obj.plotstyle, 'something went wrong with plotting!']);
                   end    
                   
                   if addcolorbar
                        h=colorbar;                           
                        h=ylabel(h, axisname); 
                        colormap jet        
                        if obj.plotconfig.latexplot
                            set(h,obj.Sfont,obj.Sint);
                        end
                   end
%                    h1=xlabel('$x_1$');
%                    h2=ylabel('$x_2$');
%                    
                   if obj.plotconfig.latexplot
%                        set(h1,obj.Sfont,obj.Sint);
%                        set(h2,obj.Sfont,obj.Sint);
%                        % axis
                       set(gca,obj.Sfont);
                   end
               end 
              % plot leaving fluxes on the boundary:
               if (any(quantity == 'vd') && obj.plotconfig.drawfluxes)
                   add_fluxes(obj_p,dat_name);
               end
               if ~obj.plotconfig.latexplot
                  title([name,title_additive]);
               end
               title_additive='';
               box on;
               %% plot additional info
                draw_geometry(obj,obj_p,dat_name);
               iplot = iplot+1;
               hold off;
           end
           if flag_drawnow
               drawnow
           end
        end          
        %%
        function draw_geometry(~,obj_p,dat_name)
           %----------------------------
           %config:
           mark_point = false;
           draw_cells = false;
           draw_connectivity = false;                %only 2d
           %connections of which points?
           AiX=1:500:obj_p.N;
           AiX = 593;
           %----------------------------
           if nargin<3
               dat_name='';
           end
           
           if obj_p.dim == 1
               % Omega
               xlim( [obj_p.Omega(1), obj_p.Omega(2)]);
               
                %%-- draw damping boundary condition
               for boun = obj_p.bc
                   a = axis;
                   if a(3) == -inf;
                       a(3)= min(obj_p.(dat_name));
                   end
                   if a(4) == inf;
                       a(4)= max(obj_p.(dat_name));
                   end
                   
                   if strcmp(boun.type,'noflow')
                        plot([boun.bp,boun.bp],a(3:4),'kx:');
                   elseif strcmp(boun.type,'nrd')
                       b = boun.damping_area;
                       rectangle('Position',[b(1), a(3), b(2)-b(1), a(4)-a(3)],...
                                 'EdgeColor',[0,1,0.5]);                   
                   elseif (strcmp (boun.type(1:3),'nrc') ...
                           || strcmp (boun.type(1:3),'nrm')...
                           || strcmp (boun.type(1:3),'nrp'))...
                           && obj_p.dim==1
                        plot([obj_p.Xj(boun.kb),obj_p.Xj(boun.kb)],a(3:4),'bx:');
                   elseif strcmp(boun.type,'cut')
                       % do nothing
                   else
                       error('No such boundary type to plot');
                   end                                         
                   axis(a);
               end
                
           elseif obj_p.dim == 2
               a = axis;
               %Omega:
               x1 = [obj_p.Omega(1,1), obj_p.Omega(2,1)];
               x2 = [obj_p.Omega(1,2), obj_p.Omega(2,1)];
               x3 = [obj_p.Omega(1,2), obj_p.Omega(2,2)];
               x4 = [obj_p.Omega(1,1), obj_p.Omega(2,2)];
               
               xlim( [obj_p.Omega(1,1), obj_p.Omega(1,2)]);
               ylim( [obj_p.Omega(2,1), obj_p.Omega(2,2)]);
               
               
               if size(a,2)>4 && ~isempty(dat_name)%->3d
%                     z =a(5);
                   z = mean(obj_p.(dat_name));
               else
                   z=0;
               end
               
               plot3( [x1(1) x2(1) x3(1) x4(1) x1(1)], [x1(2) x2(2) x3(2) x4(2) x1(2)], [z z z z z],'k' )
               
               %boundary lines:
               a = axis; %now with Omega in consideration
               for boun = obj_p.bc
                   e = [0 -1; 1 0]*boun.outer_normal'; %rotate 90degrees
                   p1= boun.bp +1000*e'; %from "minus infinity to infinity"
                   p2= boun.bp -1000*e';
                   plot3([p1(1),p2(1)],[p1(2),p2(2)],[z,z],'kx:');
               end
               axis(a);
               %% draw connectivity
               if draw_connectivity               
                   colo='gbkrmcy';
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
                   vline = [obj_p.Omega(2,1),obj_p.Omega(2,2)];
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
        %%
        function draw_points(~,obj,ind)
           hold on;
           plot(obj.Xj(ind,1),obj.Xj(ind,2),'x');            
        end
        
        
        %% save figure as eps (default)        
        function savefigure(obj,t,h,suffix)
            % t: time;
            % h: figurehandle;
            % suffix: format to save
            
            if nargin <2
               tstring='';
            else
               tstring=num2str(t);
               %remove point
               tstring = strrep(tstring, '.', '-');
            end
            if (nargin <3 || isempty(h))
                h=gcf;
            end
            if nargin <4
                suffix = obj.figure_format;
            end            
            if strcmp(suffix,'eps')
                type = 'epsc2';
            else
                type = suffix;
            end
            dir  = 'figures/';
            name = obj.plotconfig.figurename;
            if ~isempty(name)
                fullname = [dir,name,'_',tstring,'.',suffix];
%                 saveas(h,fullname,type, '-r300')    
                print(h,fullname, ['-d',type], '-r300'); %<-Save as PNG with 300 DPI
                disp (['saved figure in ',fullname]);
            else
                warning('define a figure name');
            end
        end
        
        %% %%% Write / Read Data %%% %%
        %%
        function read_hdf5(~,obj_data,filename,group)
            % obj_data: where to write the data
            % filename: what to read
            % group: which group of the hdf5 file to read (related to time)
            
            
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
            
            %Mie-Gruneisen parameter            
            if isprop(obj_data,'MG_Gammaj')
                obj_data.MG_Gammaj = readVariable('Gamma',filename,group);
                obj_data.MG_Sj = readVariable('S',filename,group);                
            end
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
          %  Gmod = readVariable('Gmod',filename,time_str);
          %  Y0 = readVariable('Y0',filename,time_str);
          %  p = readVariable('p',filename,time_str);
          %  phi = readVariable('phi',filename,time_str);
          %  tauXX = readVariable('tauXX',filename,time_str);
          %  tauXY = readVariable('tauXY',filename,time_str);
          %  tauYY = readVariable('tauYY',filename,time_str);

        end
        %%
        function write_hdf5(obj,obj_p)
            %obj_p: particle class, source of data
            
            filename = [obj.output_name,'.h5'];
            time = obj_p.t;
            function writeVariable(filename,time,name,data)
                group = ['/',num2str(time),'/',name];
                hdf5write(filename, group, data, 'WriteMode', 'append');
            end
            writeVariable(filename,time,'x',obj_p.Xj(:,1));
            writeVariable(filename,time,'u',obj_p.vj(:,1));
            writeVariable(filename,time,'c',obj_p.cj);
            writeVariable(filename,time,'e',obj_p.ej);
            writeVariable(filename,time,'c0',obj_p.c0j);
            writeVariable(filename,time,'rho',obj_p.rhoj);
            writeVariable(filename,time,'rho0',obj_p.rho0j);
            writeVariable(filename,time,'m',obj_p.mj);                        
            writeVariable(filename,time,'h',obj_p.hj);                        
            writeVariable(filename,time,'p',obj_p.pj);                        
            writeVariable(filename,time,'Gamma',obj_p.MG_Gammaj);                        
            writeVariable(filename,time,'S',obj_p.MG_Sj);  
            
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

