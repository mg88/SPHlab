%% ToDo:
% - Dissipation
% - internal energy - non-isothermal
% - rarefaction waves -> negative pressure!
% - boundary conditions (ghost particles)
% - non-reflective boundary
% - read/write hdf5 files
% - check against c++ code
% - m-scheme
% - variable h

classdef sph_particles < handle
%% SPH for HVI  - Markus Ganser - TU/e - 2015
% PARTICLES class - defines all function needed to run a sph simulation
    
%% 
    properties
        %amount of particle
        N
        dim
        %forces
        F_total
        F_int
        F_ST
        F_diss
        F_diss_art
        %density
        rho0j
        rhoj               
        rhoj_half
        drhoj               
        %pressure
        pj                  %  (p= p* = p/p0 where p0 = rho0 u0^2)
        %velocity
        vj
        vj_half
        %position
        Xj
        %connectivity ([indice of startpoint, indice of endpoint]^ammount of connections
        pij  
        %all active connection of pij in two seperate arrays 
        active_k_pij % all edges inside the cutoffradius
        Ii    %indice of point i of the ij-edge
        Ij    %indice of point j of the ij-edge
        AedgesXj;  % directed edges of each node (with respect to A_pij)
        %volume
        Vj
        %mass
        Mj
        %some indicies
        Iin
        Iall
        Iboun
        Imaterial
        %Material property 
        cj0         
        cj
        beta
        mu
            
        
        %Time
        dt
        omega %savety factor for timestepping
        tend                
        %Distances rij=xi-xj= rij * xij_h
        rij    % = bar(xij)  
        xij_h  %normalized  = hat(xij)
        nIj
        %kernel
        Wij
        dWij
        ddWij
        dWij_over_rij
        %smoothing length
        h 
        kernel    % M4 | gauss
        %Domain
        Omega
        cell_of_j
        Nc 
        Rt        
        obsolete_k_pij;
        %some flags
        distances_uptodate
        firststep
        % IO class
        IO                     
        %
        tmp
    end
 %%   
    methods
       
        %% constructor
        function obj = sph_particles(obj_scen)
            obj.IO = sph_IO(obj_scen); %initialize IO and read data if necessary
           
            obj.Omega = obj_scen.Omega;
            obj.Xj    = obj_scen.Xj;
            obj.N     = size(obj.Xj,1);
            obj.dim   = size(obj.Xj,2);

            obj.rho0j = ones(obj.N,1).*obj_scen.rho0;            
            obj.rhoj  = obj.rho0j;
            obj.cj0   = ones(obj.N,1).*obj_scen.c0;
            obj.cj    = obj.cj0;
            
            obj.beta = obj_scen.beta;
            obj.mu   = obj_scen.mu;
            
            % stimmt das?
            obj.Mj   = obj.rhoj.*obj_scen.V0particle;    
            
            %initial smoothing length
            obj.h     = obj_scen.dx * obj_scen.eta;
            obj.kernel= obj_scen.kernel;
            obj.Rt    = obj_scen.eta2* obj.h;
            
            obj.Iin   = obj_scen.Iin;
            obj.Iboun = obj_scen.Iboun;
            obj.Iall  = [obj.Iin;obj.Iboun];
            obj.Imaterial = obj_scen.Imaterial;
                                               
            obj.omega  = obj_scen.omega;
            obj.tend   = obj_scen.tend;           
            obj.vj     = obj_scen.vj;
            obj.F_int  = zeros(obj.N,obj.dim);
            obj.F_diss = zeros(obj.N,obj.dim);
            obj.F_ST   = zeros(obj.N,obj.dim);
            obj.nIj    = zeros(obj.N,obj.dim);

            obj.drhoj     = zeros(obj.N,1);
            obj.rhoj_half = zeros(obj.N,1);
            obj.vj_half   = zeros(obj.N,obj.dim);
           
            obj.pij = [];            
            obj.distances_uptodate = false;
            obj.firststep = true;
            
            obj.cell_of_j= zeros(obj.N,1);
            obj.obsolete_k_pij =[];

            %check
            checkIfInDomain(obj);
            if obj.N > 31000
                error('too many particles')
            end
            obj.tmp=0;
        end
        %%
        function start_simulation(obj)            
            disp('#### start of the simulation ####');
            %% time iteration
            t    = 0;
            icount = 1;
            tic
            ttic = tic;
            obj.IO.initialize()
            while t < obj.tend
                update_dt(obj);
                perform_timestep(obj);
                obj.IO.do(t,obj);
                t = t + obj.dt;
                %iterations per seconds
                if toc > 1 %seconds
                    disp([num2str(icount),' iter/sec ...(Ncon = ',num2str(size(obj.pij,1)),'); dt: ',num2str(obj.dt) ]);
                    icount = 1;
                    tic    
                    %disp(['tmp: ',num2str(obj.tmp)]);
                end
                icount = icount +1;     
            end
            obj.IO.finalize();
            toc(ttic)
        end
        %%        
        function perform_timestep(obj)
             update_full_step(obj)
             comp_pressure(obj);
             comp_volume(obj);
             search_neighbours(obj);
             update_kernel(obj)
             comp_forces(obj)
             comp_dRho(obj)
             update_half_step(obj)
             update_postion(obj)           
        end     
        %%
        function update_dt(obj)  %absolute value ?!
            if obj.dim == 1
                obj.dt = obj.omega * min(obj.h./(obj.cj + obj.vj));
            else
                obj.dt = obj.omega * min(obj.h./(obj.cj + sum(obj.vj(:,1).^2 + obj.vj(:,2).^2,2).^(0.5)));
            end
        end
        
        %% neighbour search -  main 
        function search_neighbours(obj)
           if obj.firststep
               initial_search_neighbours(obj);
           else
               update_neighbours(obj);
               %initial_search_neighbours(obj);
           end           
        end        
        %%
        function initial_search_neighbours(obj) %upper right and lower left cells are not included
            
            %% create cell structure
            obj.Nc = floor((obj.Omega(:,2)-obj.Omega(:,1))/obj.Rt);
            if obj.dim ==1
                obj.cell_of_j = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                cellshift  = 1;
                Ncell_search = obj.Nc-1;
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                obj.cell_of_j = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(1));    %convert to 1d counting
                cellshift = [1,obj.Nc(1)-1,obj.Nc(1),obj.Nc(1)+1]; %upper-left side
                Ncell_search = (obj.Nc(1)*obj.Nc(2) - (obj.Nc(1)+1));
            end
            %obj.Kcells=sparse(c1d,1:obj.N,1);              %matrix includes a 1-entry if a particle (colum) is in a cell (row)            
            
            if any(obj.cell_of_j > Ncell_search)
                %keyboard
                error(' - particles not in cell structure! - ')
            end
            
            %% create connectivity list:
            obj.pij=[];
            kp=1;
            for kcell = 1:Ncell_search %loop over all cells except the upper row                
                j_in_Cell = find(obj.cell_of_j==kcell);
                if size(j_in_Cell,1)>0   
                   i_in_SurroundingCell = find(ismember(obj.cell_of_j,kcell+cellshift)); 
                   % i_in_SurroundingCell =
                   % find(ismembc(obj.cell_of_j,kcell+cellshift));
                   % %undocumented version which might be faster (second
                   % argument has to be sorted
                    i_all = [j_in_Cell;i_in_SurroundingCell];
                    nParticle_all = length(i_all);
                    for kj = 1:(length(j_in_Cell)) %loop over particles in the cell                    
                        %% with cutoff:
    %                     allEdges= ones(nParticle_all-ki,1)*obj.Xj(i_all(ki),:) - obj.Xj(i_all(ki+1:end),:);
    %                     d=sum(allEdges.^2,2);
    %                     [Y,I] = sort(abs(d)); %sort
    %                     neighbors_of_j = I(Y<Rtsquare); 
                        %% without cutoff
                       neighbors_of_j = ((kj+1):nParticle_all)'-kj; 
                        %%
                        NN=size(neighbors_of_j,1);
                        obj.pij(kp:kp+NN-1,:)=...
                           [j_in_Cell(kj)*ones(NN,1),i_all(neighbors_of_j+kj)];
                        kp=kp+NN; 
                    end
                end    
            end
            %% create adjacency matrix nodes <-> edges
            obj.AedgesXj = sparse(obj.N,size(obj.pij,1));
            for kk = 1:obj.N
                con = sparse(obj.pij==kk);
                obj.AedgesXj(kk,con(:,1))= 1;  %outgoing edges
                obj.AedgesXj(kk,con(:,2))= -1; % ingoing edges                  
            end 
        end
        %%
        function Ashift = lockup_cellshift(obj,jshift)
           if obj.dim==1
              switch jshift
                  case 1
                      Ashift = 1;
                  case -1
                      Ashift = -1;
                  otherwise
                        error('A particel skipped a cell!')
              end
                     
           else
               switch jshift
                   case 1 %right
                        Ashift =  1+[obj.Nc(1),0,-obj.Nc(1)];
                   case -1 %left
                        Ashift = -1+[obj.Nc(1),0,-obj.Nc(1)];
                   case obj.Nc(1)%up
                        Ashift = obj.Nc(1)+[1,0,-1];
                   case obj.Nc(1)+1 % upper right
                        Ashift = [-obj.Nc(1)+1,1,1+obj.Nc(1),obj.Nc(1),obj.Nc(1)-1];
                   case obj.Nc(1)-1 % upper left
                       Ashift =  [-obj.Nc(1)-1,-1,-1+obj.Nc(1),obj.Nc(1),obj.Nc(1)+1];
                   case -obj.Nc(1) %down
                        Ashift = -obj.Nc(1)+[1,0,-1];
                   case -obj.Nc(1)-1 % down left
                        Ashift = [-obj.Nc(1)-1,-obj.Nc(1),-obj.Nc(1)+1,-1,obj.Nc(1)-1];
                   case -obj.Nc(1)+1 % down right
                        Ashift = [-obj.Nc(1)-1,-obj.Nc(1),-obj.Nc(1)+1,+1,obj.Nc(1)+1];
                   otherwise
                        error('A particel skipped a cell!')
               end
           end   
        end     
        %%
        function update_neighbours (obj) 
            %% create cell structure
            if obj.dim == 1
                cell_of_j_new = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle      
                if any(cell_of_j_new > obj.Nc)  || any(cell_of_j_new < 0)
                    %keyboard
                    error(' - particles not in cell structure anymore! - ')
                end
           
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                cell_of_j_new = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(1));    %convert to 1d counting
                if any(cell_of_j_new > obj.Nc(1)*obj.Nc(2))  || any(cell_of_j_new < 0)
                    %keyboard
                    error(' - particles not in cell structure anymore! - ')
                end
            end
              
                        

            
            cell_shift_of_j = cell_of_j_new - obj.cell_of_j;
            j_changes_cell = find(cell_shift_of_j); %alternative: a-b ~=0
                  
            %size(j_changes_cell,1)
            
            %% loop over all particels which have changed the cell
            for j = j_changes_cell'  
                cellshift = lockup_cellshift(obj,cell_shift_of_j(j)); %cells ahead in moving direction
                
                %find obsolete points (all particels in cells which are no
                %neighbours any more - right now, one behind ...|X|j|o|o|...
                cells_to_look_at = obj.cell_of_j(j) - cellshift;
                i_in_obs_neighb_cell = sparse(double(ismember(obj.cell_of_j,cells_to_look_at)));
                
                %find the belonging connectivities
                i_in_obs_neighb_cell(j,1) = 5;  % just a number,to obtain an unsymmetric connection - abs(+-5 -+1)=4
                temp = (i_in_obs_neighb_cell'*obj.AedgesXj);
                pij_obs = (abs(temp)==4);
                obj.obsolete_k_pij = [obj.obsolete_k_pij;...
                                         find(pij_obs')];
                    
                % move j to new cell (each particel at one time for
                % consistency)     % ...|o|j|o|o|... ->  ...|o|o|j|o|...
                obj.cell_of_j(j) = cell_of_j_new(j);

                %find new points in cells ahead ...|o|o|j|X|...
                cells_to_look_at = obj.cell_of_j(j) + cellshift;
                i_in_new_neighb_cell = find(ismember(obj.cell_of_j,cells_to_look_at));
                nCon_new = size(i_in_new_neighb_cell,1);                            
                                                                
                % clear all obsolete connectivities of node j in the
                % adjacency matrix nodes <-> edges
                obj.AedgesXj(:,obj.obsolete_k_pij) = 0;
                % http://de.mathworks.com/matlabcentral/newsreader/view_thread/283923
                             
                if nCon_new > 0
                    %%new connectivies               
                    pij_new = [j*ones(nCon_new,1),i_in_new_neighb_cell]; 
                    nCon_obsolete = length(obj.obsolete_k_pij);
                    if nCon_new <= nCon_obsolete
                         k_pij_new =  obj.obsolete_k_pij(1:nCon_new);
                         obj.pij(k_pij_new,:) = pij_new;   
                         obj.obsolete_k_pij = obj.obsolete_k_pij(nCon_new+1:end);    
                         obj.pij(obj.obsolete_k_pij,:) =...;
                            ones(nCon_obsolete-nCon_new,1)*[j,j];  %set a self connection 
                    else
                         %write on present connectivities
                         obj.pij(obj.obsolete_k_pij,:)= pij_new (1:nCon_obsolete,:);
                         %and create new ones
                         n_pij_new = nCon_new - nCon_obsolete;
                         obj.pij = [obj.pij;
                                     pij_new(nCon_obsolete+1:nCon_new,:)];
                         n_Apij    = size(obj.pij,1);
                         k_pij_new = [obj.obsolete_k_pij;...
                                     ((n_Apij-n_pij_new+1):n_Apij)'];
                         % enlarge adjacency matrix 
                         obj.AedgesXj(obj.N,n_Apij)=0;
                    end

                    %% update adjacency matrix 
                    %second, write new connectivities in adjacency matrix
                    obj.AedgesXj (j , k_pij_new) = 1; %outgoing
                    obj.AedgesXj (i_in_new_neighb_cell,k_pij_new) = -1*speye(nCon_new); %ingoing  
                end
            end              
        end                        
        %%
        function comp_distances(obj)
            obj.xij_h=(obj.Xj(obj.pij(:,1),:)-obj.Xj(obj.pij(:,2),:));
            obj.rij= sum(obj.xij_h.^2,2).^0.5;
            obj.xij_h = obj.xij_h./(obj.rij*ones(1,obj.dim));  
            
            % flag particle which are within the cutoff radius           ^       
            temp = (obj.rij < obj.Rt).*(obj.rij > 0);
            obj.active_k_pij = logical(temp);
            obj.Ii = obj.pij(logical(temp),1);
            obj.Ij = obj.pij(logical(temp),2);
            
            obj.distances_uptodate = true;
        end  
        %%                 
        function comp_normal(obj)  %violeau422    
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            %compute the flux
             qA_nIij=zeros(size(obj.pij,1),obj.dim);
             qA_nIij(obj.active_k_pij,:) = (obj.Vj(obj.Ij,:).*obj.dWij(:)... Vj*dWij*hat(rij)
                 *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:);
             %spread fluxes to the nodes
             for dimension=1:obj.dim
                obj.nIj(:,dimension)  = obj.AedgesXj * qA_nIij(:,dimension);
             end
             %normalize:           
             norma=(sum(obj.nIj.^2,2)).^0.5;
             I=abs(norma)<eps;
             obj.nIj=obj.nIj./(norma*ones(1,obj.dim));
             obj.nIj(I,:)=0;    %take care for the almost zero values - set them zero                         
        end
         
        %% - choose a kernel function
        function update_kernel(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            if strcmp(obj.kernel,'gauss')
                update_kernel_gauss(obj);
            elseif strcmp(obj.kernel,'M4')
                update_kernel_M4(obj);  
            else
                error([obj.kernel,' as kernel is not implemented yet']);
            end
        end
        %% Gauss kernel
        function update_kernel_gauss(obj)
            sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
            obj.Wij   = 1/(obj.h^obj.dim)* sigma(obj.dim) *...
                        exp(-(obj.rij(obj.active_k_pij)/obj.h).^2);
            obj.dWij  = obj.Wij.* -2.*obj.rij(obj.active_k_pij)/obj.h^2;
            obj.dWij_over_rij = obj.Wij.* -2./obj.h^2;
            obj.ddWij = -2/obj.h^2*(obj.dWij.* obj.rij(obj.active_k_pij) + obj.Wij );
        end
        %% M4 spline kernel
        function update_kernel_M4(obj)
            sigma = [2/3; 10/(7*pi); 1/pi] ./ (obj.h.^obj.dim); %normalisation constant
            r = obj.rij(obj.active_k_pij) / obj.h;
     %       r = linspace(0,2,100);
            obj.Wij   = sigma (obj.dim) .*...
                (r<2)  .*(...
                1/4*(2-r).^3   ...
                - (1-r).^3 .* (r<1)...
                );
            
            obj.dWij   = 1/obj.h * sigma (obj.dim) .*...
                (r<2)  .*(...
                -3/4*(2-r).^2  ...
                + 3*(1-r).^2 .* (r<1)...
                );
                        
            obj.dWij_over_rij = obj.dWij ./ obj.rij(obj.active_k_pij);
            
            obj.ddWij = 1/obj.h^2 * sigma (obj.dim) .* ...
                (r<2)  .*(...
                +6/4*(2-r)  ...
                - 6*(1-r) .* (r<1)...
                ); 
    %       plot(r,obj.A_Wij,r,obj.A_dWij, r,obj.A_ddWij);
    %       keyboard
        end              
        %%
        function update_full_step(obj) 
            if ~obj.firststep
                obj.rhoj(obj.Iall) = obj.rhoj_half(obj.Iall)+ 0.5*obj.dt * obj.drhoj(obj.Iall);
                obj.vj(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + 0.5*obj.dt *...
                    (obj.F_int(obj.Iin,:)+obj.F_ST(obj.Iin,:)+obj.F_diss(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
            end
        end
        %%
        function update_half_step(obj) 
            if obj.firststep
                 obj.rhoj_half(obj.Iall) = obj.rhoj(obj.Iall)  + 0.5*obj.dt*obj.drhoj(obj.Iall);
                 obj.vj_half(obj.Iin,:)  = obj.vj(obj.Iin,:) + 0.5*obj.dt * ...
                    (obj.F_total(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
                obj.firststep=false;
            else
                obj.rhoj_half(obj.Iall) = obj.rhoj_half(obj.Iall)+ obj.dt * obj.drhoj(obj.Iall);
                obj.vj_half(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + obj.dt *...
                    (obj.F_total(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
                %+ones(size(Iin))*g_ext);
            end
        end
        %%
        function update_postion(obj) 
             obj.Xj(obj.Iin,:)      = obj.Xj(obj.Iin,:) + obj.dt * obj.vj_half(obj.Iin,:);
             obj.distances_uptodate = false;
        end

        %%
        function comp_volume(obj)
            obj.Vj = obj.Mj(:)./obj.rhoj(:);
        end        
        %%
        function comp_pressure(obj) %switch
            comp_pressure_isothermal(obj)
        end
        
        function comp_pressure_isothermal(obj)
            obj.pj = obj.cj0 .* obj.cj0 .*(obj.rhoj - obj.rho0j);
        end
        %%
        function comp_dRho(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            obj.drhoj  = zeros(size(obj.drhoj)); 
            qrho_ij = zeros(size(obj.pij,1),1);
            qrho_ij(obj.active_k_pij,:) = obj.rhoj(obj.Ii,:).*obj.Vj(obj.Ij,:).*... %density flux
                sum(... %scalar product v*n 
                (obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)) .*...
                (obj.dWij*ones(1,obj.dim)).*obj.xij_h((obj.active_k_pij),:)...
                ,2);
            %add up all the corresponding density flux in each node
            obj.drhoj = abs(obj.AedgesXj) * qrho_ij;      
        end        
        %%
        function comp_forces(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            comp_Fint(obj);
            comp_Fdiss_art(obj)
            if obj.mu ~=0
                comp_Fdiss(obj);
            end
            if obj.beta ~= 0
                comp_F_ST(obj);
            end
            obj.F_total = obj.F_int + obj.F_diss + obj.F_diss_art + obj.F_ST;
        end
        %%
        function comp_Fint(obj) %internal forces  
            %compute the flux
            qF_int_ij = zeros(size(obj.pij,1),obj.dim);
            qF_int_ij(obj.active_k_pij,:)   = - ((obj.Vj(obj.Ii,:) .* obj.Vj(obj.Ij,:) .*... %1=i; 2=j
                                 (obj.pj(obj.Ii,:)+obj.pj(obj.Ij,:)) .* ... 
                                  obj.dWij)...
                              *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:);
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_int(:,dimension)  = obj.AedgesXj * qF_int_ij(:,dimension);
            end
%             if any(isnan(obj.F_int))
%                keyboard
%             end
        end
        %%
        function comp_Fdiss(obj)      %page415-violeau | particle-friction   
            %compute the flux
            qF_diss_ij = zeros(size(obj.pij,1),obj.dim);
            qF_diss_ij(obj.active_k_pij,:)  = 2*(obj.dim+2)*obj.mu*((...
                obj.Vj(obj.Ii,:) .* obj.Vj(obj.Ij,:).*... Vi*Vj
                         sum(...
                         ((obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)).*... $(vi-vj)*hat(xij)
                         obj.xij_h(obj.active_k_pij,:))...
                         ,2).*...
                         obj.dWij)...  %dWij
                         *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:); %hat(xij)
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_diss(:,dimension)  = obj.AedgesXj * qF_diss_ij(:,dimension);
            end
        end
        %%
        function comp_Fdiss_art(obj)      %dissipative velocity - Iason %h constant!
            %compute the flux
            qF_diss_art_ij = zeros(size(obj.pij,1),obj.dim);
            alpha_diss = 1;
            beta_diss  = 2;
            vijxijh = sum(...
                   ((obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:))).*... (vi-vj)*hat(xij)
                     obj.xij_h(obj.active_k_pij,:)...
                     ,2);
             v_sig = (0.5 *(obj.cj(obj.Ii)+obj.cj(obj.Ij))...
                   - 0.5*beta_diss*vijxijh)...
                   .*(vijxijh<=0); 
               
             qF_diss_art_ij(obj.active_k_pij,:)  = ...
                 ((obj.Mj(obj.Ii) .* obj.Mj(obj.Ij).* ...
                 alpha_diss .*v_sig .* vijxijh) ./...
                 (0.5 *(obj.rhoj(obj.Ii) + obj.rhoj(obj.Ij)))... bar(rho)
                 .*obj.dWij)...   %(h constant!)
                 *ones(1,obj.dim)...
                 .*obj.xij_h(obj.active_k_pij,:);
             
            for dimension=1:obj.dim
                obj.F_diss_art(:,dimension)  = obj.AedgesXj * qF_diss_art_ij(:,dimension);
            end
        end        
        %%
        function comp_F_ST(obj)  %Violeau423       %minus vor beta?       
             comp_normal(obj)
            %compute the flux
             qF_ST_ij = zeros(size(obj.pij,1),obj.dim);
             qF_ST_ij(obj.active_k_pij,:)  =  - obj.beta* ((obj.Vj(obj.Ii,:) .* obj.Vj(obj.Ij,:)) ... Vi*Vj
                     *ones(1,obj.dim)).*...
                    (((obj.ddWij(:) - obj.dWij_over_rij(:)).*...  ddWij-dWij/rij
                    sum(...
                    (obj.xij_h(obj.active_k_pij,:).*(obj.nIj(obj.Ii,:)-obj.nIj(obj.Ij,:)))... xij_h*(ni-nj)
                    ,2)...
                    *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:)... xij_h
                    +...
                     (obj.dWij_over_rij(:)*ones(1,obj.dim)).*... %dWij/rij
                     (obj.nIj(obj.Ii,:)-obj.nIj(obj.Ij,:))); % ni-nj
            
            
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_ST(:,dimension)  = obj.AedgesXj * qF_ST_ij(:,dimension);
            end            
        end
        %%
        function checkIfInDomain(obj)
            if (any(any(obj.Xj-  ones(obj.N,1)* obj.Omega(:,1)'<0)) || any(any( ones(obj.N,1)*obj.Omega(:,2)' - obj.Xj <0)))
                error('some points are out of the domain');
            end
        end              
    end
end
