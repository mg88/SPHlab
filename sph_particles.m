%% ToDo:
% - Dissipation
% - internal energy - non-isothermal
% - rarefaction waves!
% - boundary conditions (ghost particles)
% - non-reflective boundary
% - read/write hdf5 files
% - check against c++ code

classdef sph_particles < handle
%% SPH for HVI  - Markus Ganser - TU/e - 2015
% PARTICLES class - defines all function needed to run a sph simulation
    
%% 
    properties
        %amount of particle
        N
        dim
        %forces
        F_int
        F_ST
        F_diss
        %density
        Rhoj
        Rhoj_half
        dRhoj
        %pressure
        pj
        %velocity
        vj
        vj_half
        %position
        Xj
        %connectivity
        A_pij
        A_edgesXj;  % directed edges of each node (with respect to A_pij)
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
        beta
        mu
        Ca       
        %Time
        dt
        tend                
        %Distances
        A_rij
        A_eij
        A_nIj
        %kernel
        A_Wij
        A_dWij
        A_ddWij
        A_dWij_over_rij
        %smoothing length
        h 
        %Domain
        Omega
        cj %cell division structure
        cell_of_j
        Ncells
        Nc  %both neessary?
        Rt
        active_k_pij
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
        %function obj = particles(obj_geo,rho0,h,Iin,Iboun)
        function obj = sph_particles(obj_scen)
            obj.IO = sph_IO(obj_scen);
            if obj_scen.read_file
               obj_scen = obj.IO.read(obj_scen.input_name);
            end
            
            
            obj.Omega = obj_scen.Omega;
            obj.Xj    = obj_scen.obj_geo.Xj;
            obj.N     = size(obj.Xj,1);
            obj.dim   = size(obj.Xj,2);

            obj.Rhoj = ones(obj.N,1).*obj_scen.rho0;
            obj.Mj   = obj.Rhoj.*obj_scen.obj_geo.V0particle;    
            
            %initial smoothing length
            obj.h     = obj_scen.dx * obj_scen.eta;
            obj.Rt    = obj_scen.eta2* obj.h;
            
            obj.Iin   = obj_scen.Iin;
            obj.Iboun = obj_scen.Iboun;
            obj.Iall  = [obj.Iin;obj.Iboun];
            obj.Imaterial = obj_scen.Imaterial;
            
            obj.beta = obj_scen.beta;
            obj.mu   = obj_scen.mu;
            obj.Ca   = obj_scen.Ca;
            
            obj.dt   = obj_scen.dt;
            obj.tend = obj_scen.tend;
            
            obj.vj     = obj_scen.obj_geo.vj;
            obj.F_int  = zeros(obj.N,obj.dim);
            obj.F_diss = zeros(obj.N,obj.dim);
            obj.F_ST   = zeros(obj.N,obj.dim);
            obj.A_nIj  = zeros(obj.N,obj.dim);

            obj.dRhoj     = zeros(obj.N,1);
            obj.Rhoj_half = zeros(obj.N,1);
            obj.vj_half   = zeros(obj.N,obj.dim);
           
            obj.A_pij = [];
            obj.distances_uptodate = false;
            obj.firststep = true;
            
            obj.cell_of_j= zeros(obj.N,1);
            obj.obsolete_k_pij =[];

            %check
            checkIfInDomain(obj);
            if obj.N > 30000
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
            ttic = tic;
            obj.IO.initialize()
            while t < obj.tend
                perform_timestep(obj);
                obj.IO.do(t,obj);
                t = t + obj.dt;
                
                %iterations per seconds
                if toc > 1 %seconds
                    disp([num2str(icount),' iter/sec ...(Ncon = ',num2str(size(obj.A_pij,1)),')']);
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
            obj.Nc = floor(obj.Omega/obj.Rt);
            if obj.dim ==1
                obj.cell_of_j = floor((ones(obj.N,1)*obj.Nc./obj.Omega) .* obj.Xj)+1; %in which sector is the particle
                obj.Ncells = obj.Nc;
                cellshift  = 1;
                Ncell_search = obj.Ncells-1;
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc./obj.Omega) .* obj.Xj)+1; %in which sector is the particle
                obj.cell_of_j = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(2));    %convert to 1d counting
                obj.Ncells=obj.Nc(1)*obj.Nc(2);
                cellshift = [1,obj.Nc(2)-1,obj.Nc(2),obj.Nc(2)+1]; %upper-left side
                Ncell_search = (obj.Ncells - (obj.Nc(2)+1));
            end
            %obj.Kcells=sparse(c1d,1:obj.N,1);              %matrix includes a 1-entry if a particle (colum) is in a cell (row)            
            
            %% create connectivity list:
            obj.A_pij=[];
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
                        obj.A_pij(kp:kp+NN-1,:)=...
                           [j_in_Cell(kj)*ones(NN,1),i_all(neighbors_of_j+kj)];
                        kp=kp+NN; 
                    end
                end    
            end
            
            %% create adjacency matrix nodes <-> edges
            obj.A_edgesXj = sparse(obj.N,size(obj.A_pij,1));
            for kk = 1:obj.N
                con = sparse(obj.A_pij==kk);
                obj.A_edgesXj(kk,con(:,1))= 1;  %outgoing edges
                obj.A_edgesXj(kk,con(:,2))= -1; % ingoing edges                  
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
                        Ashift =  1+[obj.Nc,0,-obj.Nc(2)];
                   case -1 %left
                        Ashift = -1+[obj.Nc,0,-obj.Nc(2)];
                   case obj.Nc(2)%up
                        Ashift = obj.Nc(2)+[1,0,-1];
                   case obj.Nc(2)+1 % upper right
                        Ashift = [-obj.Nc(2)+1,1,1+obj.Nc(2),obj.Nc(2),obj.Nc(2)-1];
                   case obj.Nc(2)-1 % upper left
                       Ashift = [-obj.Nc(2)-1,-1,-1+obj.Nc(2),obj.Nc(2),obj.Nc(2)+1];
                   case -obj.Nc(2) %down
                        Ashift = -obj.Nc(2)+[1,0,-1];
                   case -obj.Nc(2)-1 % down left
                        Ashift = [-obj.Nc(2)-1,-obj.Nc(2),-obj.Nc(2)+1,-1,obj.Nc(2)-1];
                   case -obj.Nc(2)+1 % down right
                        Ashift = [-obj.Nc(2)-1,-obj.Nc(2),-obj.Nc(2)+1,+1,obj.Nc(2)+1];
                   otherwise
                        error('A particel skipped a cell!')
               end
           end   
        end     
        %%
        function update_neighbours (obj) 
            %% create cell structure
            if obj.dim ==1
                cell_of_j_new = floor((ones(obj.N,1)*obj.Nc./obj.Omega) .* obj.Xj)+1; %in which sector is the particle
                obj.Ncells=obj.Nc;
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc./obj.Omega) .* obj.Xj)+1; %in which sector is the particle
                cell_of_j_new = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(2));    %convert to 1d counting
                obj.Ncells=obj.Nc(1)*obj.Nc(2);
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
                temp = (i_in_obs_neighb_cell'*obj.A_edgesXj);
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
                obj.A_edgesXj(:,obj.obsolete_k_pij) = 0;
                % http://de.mathworks.com/matlabcentral/newsreader/view_thread/283923
                             
                if nCon_new > 0
                    %%new connectivies               
                    pij_new = [j*ones(nCon_new,1),i_in_new_neighb_cell]; 
                    nCon_obsolete = length(obj.obsolete_k_pij);
                    if nCon_new <= nCon_obsolete
                         k_pij_new =  obj.obsolete_k_pij(1:nCon_new);
                         obj.A_pij(k_pij_new,:) = pij_new;   
                         obj.obsolete_k_pij = obj.obsolete_k_pij(nCon_new+1:end);    
                         obj.A_pij(obj.obsolete_k_pij,:) =...;
                            ones(nCon_obsolete-nCon_new,1)*[j,j];  %set a self connection 
                    else
                         %write on present connectivities
                         obj.A_pij(obj.obsolete_k_pij,:)= pij_new (1:nCon_obsolete,:);
                         %and create new ones
                         n_pij_new = nCon_new - nCon_obsolete;
                         obj.A_pij = [obj.A_pij;
                                     pij_new(nCon_obsolete+1:nCon_new,:)];
                         n_Apij    = size(obj.A_pij,1);
                         k_pij_new = [obj.obsolete_k_pij;...
                                     ((n_Apij-n_pij_new+1):n_Apij)'];
                         % enlarge adjacency matrix 
                         obj.A_edgesXj(obj.N,n_Apij)=0;
                    end

                    %% update adjacency matrix 
                    %second, write new connectivities in adjacency matrix
                    obj.A_edgesXj (j , k_pij_new) = 1; %outgoing
                    obj.A_edgesXj (i_in_new_neighb_cell,k_pij_new) = -1*speye(nCon_new); %ingoing  
                end
            end              
        end                
        %%                 
        function comp_normal(obj)  %violeau422    
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            %compute the flux
             qA_nIij=zeros(size(obj.A_pij,1),obj.dim);
             qA_nIij(obj.active_k_pij,:) = (obj.Vj(obj.A_pij(obj.active_k_pij,2),:).*obj.A_dWij(:)... Vj*dWij*eij
                 *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:);
             %spread fluxes to the nodes
             for dimension=1:obj.dim
                obj.A_nIj(:,dimension)  = obj.A_edgesXj * qA_nIij(:,dimension);
             end
             %normalize:           
             norma=(sum(obj.A_nIj.^2,2)).^0.5;
             I=abs(norma)<eps;
             obj.A_nIj=obj.A_nIj./(norma*ones(1,obj.dim));
             obj.A_nIj(I,:)=0;    %take care for the almost zero values - set them zero                         
        end
        %%
        function comp_distances(obj)
            obj.A_eij=(obj.Xj(obj.A_pij(:,1),:)-obj.Xj(obj.A_pij(:,2),:));
            obj.A_rij= sum(obj.A_eij.^2,2).^0.5;
            obj.A_eij = obj.A_eij./(obj.A_rij*ones(1,obj.dim));           
            % flag particle which are within the cutoff radius           ^       
            temp = logical(obj.A_rij < obj.Rt).*(obj.A_rij > 0);
            obj.active_k_pij = logical(temp);
            
            obj.distances_uptodate = true;
        end           
        %% - choose a kernel function
        function update_kernel(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
             %update_kernel_gauss(obj);
             update_kernel_M4(obj);  
        end
        %% Gauss kernel
        function update_kernel_gauss(obj)
            sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
            obj.A_Wij   = 1/(obj.h^obj.dim)* sigma(obj.dim) *...
                        exp(-(obj.A_rij(obj.active_k_pij)/obj.h).^2);
            obj.A_dWij  = obj.A_Wij.* -2.*obj.A_rij(obj.active_k_pij)/obj.h^2;
            obj.A_dWij_over_rij = obj.A_Wij.* -2./obj.h^2;
            obj.A_ddWij = -2/obj.h^2*(obj.A_dWij.* obj.A_rij(obj.active_k_pij) + obj.A_Wij );
        end
        %% M4 spline kernel
        function update_kernel_M4(obj)
            sigma = [2/3; 10/(7*pi); 1/pi] ./ (obj.h.^obj.dim); %normalisation constant
            r = obj.A_rij(obj.active_k_pij) / obj.h;
     %       r = linspace(0,2,100);
            obj.A_Wij   = sigma (obj.dim) .*...
                (r<2)  .*(...
                1/4*(2-r).^3   ...
                - (1-r).^3 .* (r<1)...
                );
            
            obj.A_dWij   = 1/obj.h * sigma (obj.dim) .*...
                (r<2)  .*(...
                -3/4*(2-r).^2  ...
                + 3*(1-r).^2 .* (r<1)...
                );
                        
            obj.A_dWij_over_rij = obj.A_dWij ./ obj.A_rij(obj.active_k_pij);
            
            obj.A_ddWij = 1/obj.h^2 * sigma (obj.dim) .* ...
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
                obj.Rhoj(obj.Iall) = obj.Rhoj_half(obj.Iall)+ 0.5*obj.dt * obj.dRhoj(obj.Iall);
                obj.vj(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + 0.5*obj.dt *...
                    (obj.F_int(obj.Iin,:)+obj.F_ST(obj.Iin,:)+obj.F_diss(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
            end
        end
        %%
        function update_half_step(obj) 
            if obj.firststep
                 obj.Rhoj_half(obj.Iall) = obj.Rhoj(obj.Iall)  + 0.5*obj.dt*obj.dRhoj(obj.Iall);
                 obj.vj_half(obj.Iin,:)  = obj.vj(obj.Iin,:) + 0.5*obj.dt * ...
                    (obj.F_int(obj.Iin,:)+obj.F_ST(obj.Iin,:)+obj.F_diss(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
                obj.firststep=false;
            else
                obj.Rhoj_half(obj.Iall) = obj.Rhoj_half(obj.Iall)+ obj.dt * obj.dRhoj(obj.Iall);
                obj.vj_half(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + obj.dt *...
                    (obj.F_int(obj.Iin,:)+obj.F_ST(obj.Iin,:)+obj.F_diss(obj.Iin,:))./(obj.Mj(obj.Iin)*ones(1,obj.dim));
                %+ones(size(Iin))*g_ext);
            end
        end
        %%
        function update_postion(obj) 
             obj.Xj(obj.Iin,:)      = obj.Xj(obj.Iin,:) + obj.dt * obj.vj_half(obj.Iin,:);
             obj.distances_uptodate = false;
        end
        %%
        function comp_dRho(obj) %evtl. reshape instead of ones(1,dim)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            obj.dRhoj  = zeros(size(obj.dRhoj)); 
            qrho_ij = zeros(size(obj.A_pij,1),1);
            qrho_ij(obj.active_k_pij,:) = obj.Rhoj(obj.A_pij((obj.active_k_pij),1),:).*obj.Vj(obj.A_pij((obj.active_k_pij),2),:).*... %density flux
                sum(... %scalar product v*n 
                (obj.vj(obj.A_pij((obj.active_k_pij),1),:)-obj.vj(obj.A_pij((obj.active_k_pij),2),:)) .*...
                (obj.A_dWij*ones(1,obj.dim)).*obj.A_eij((obj.active_k_pij),:)...
                ,2);
            %add up all the corresponding density flux in each node
            obj.dRhoj = abs(obj.A_edgesXj) * qrho_ij;      
        end
        %%
        function comp_volume(obj)
            obj.Vj=obj.Mj(:)./obj.Rhoj(:);
        end        
        %%
        function comp_pressure(obj)
            obj.pj = 1./obj.Ca.*(obj.Rhoj-1);
        end
        %%
        function comp_forces(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            comp_Fint(obj);
            if obj.mu ~=0
                comp_Fdiss(obj);
            end
            if obj.beta ~= 0
                comp_F_ST(obj);
            end
        end
        %%
        function comp_Fint(obj) %internal forces  
            %compute the flux
            qF_int_ij = zeros(size(obj.A_pij,1),obj.dim);
            qF_int_ij(obj.active_k_pij,:)   = - ((obj.Vj(obj.A_pij(obj.active_k_pij,1),:) .* obj.Vj(obj.A_pij(obj.active_k_pij,2),:) .*... %1=i; 2=j
                                 (obj.pj(obj.A_pij(obj.active_k_pij,1),:)+obj.pj(obj.A_pij(obj.active_k_pij,2),:)) .* ... 
                                  obj.A_dWij)...
                              *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:);
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_int(:,dimension)  = obj.A_edgesXj * qF_int_ij(:,dimension);
            end
            if any(isnan(obj.F_int))
               keyboard
            end

        end
        %%
        function comp_Fdiss(obj)      %page415-violeau | particle-friction   
            %compute the flux
            qF_diss_ij = zeros(size(obj.A_pij,1),obj.dim);
            qF_diss_ij(obj.active_k_pij,:)  = 2*(obj.dim+2)*obj.mu*((...
                obj.Vj(obj.A_pij(obj.active_k_pij,1),:) .* obj.Vj(obj.A_pij(obj.active_k_pij,2),:).*... Vi*Vj
                         sum(...
                         ((obj.vj(obj.A_pij(obj.active_k_pij,1),:)-obj.vj(obj.A_pij(obj.active_k_pij,2),:)).*... $(vi-vj)*eij
                         obj.A_eij(obj.active_k_pij,:))...
                         ,2).*...
                         obj.A_dWij)...  %dWij
                         *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:); %eij
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_diss(:,dimension)  = obj.A_edgesXj * qF_diss_ij(:,dimension);
            end
        end
        %%
        function comp_Fdiss_art(obj)      %dissipative velocity - Iason
            %compute the flux
            qF_diss_ij = zeros(size(obj.A_pij,1),obj.dim);
           % alpha = 1;
           % beta  = 2;
            
%             v_sig = -0.5*beta((obj.vj(obj.A_pij(obj.active_k_pij,1),:)-(obj.A_pij(obj.active_k_pij,2),:)).*...
%                     *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:); %eij
            
            qF_diss_ij(obj.active_k_pij,:)  = 2*(obj.dim+2)*obj.mu*((...
                obj.Vj(obj.A_pij(obj.active_k_pij,1),:) .* obj.Vj(obj.A_pij(obj.active_k_pij,2),:).*... Vi*Vj
                         sum(...
                         ((obj.vj(obj.A_pij(obj.active_k_pij,1),:)-obj.vj(obj.A_pij(obj.active_k_pij,2),:)).*... $(vi-vj)*eij
                         obj.A_eij(obj.active_k_pij,:))...
                         ,2).*...
                         obj.A_dWij)...  %dWij
                         *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:); %eij
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_diss(:,dimension)  = obj.A_edgesXj * qF_diss_ij(:,dimension);
            end
        end        
        %%
        function comp_F_ST(obj)  %Violeau423       %minus vor beta?       
             comp_normal(obj)
            %compute the flux
             qF_ST_ij = zeros(size(obj.A_pij,1),obj.dim);
             qF_ST_ij(obj.active_k_pij,:)  =  - obj.beta* ((obj.Vj(obj.A_pij(obj.active_k_pij,1),:) .* obj.Vj(obj.A_pij(obj.active_k_pij,2),:)) ... Vi*Vj
                     *ones(1,obj.dim)).*...
                    (((obj.A_ddWij(:) - obj.A_dWij_over_rij(:)).*...  ddWij-dWij/rij
                    sum(...
                    (obj.A_eij(obj.active_k_pij,:).*(obj.A_nIj(obj.A_pij(obj.active_k_pij,1),:)-obj.A_nIj(obj.A_pij(obj.active_k_pij,2),:)))... eij*(ni-nj)
                    ,2)...
                    *ones(1,obj.dim)).*obj.A_eij(obj.active_k_pij,:)... eij
                    +...
                     (obj.A_dWij_over_rij(:)*ones(1,obj.dim)).*... %dWij/rij
                     (obj.A_nIj(obj.A_pij(obj.active_k_pij,1),:)-obj.A_nIj(obj.A_pij(obj.active_k_pij,2),:))); % ni-nj
            
            
            %spread fluxes to the nodes
            for dimension=1:obj.dim
                obj.F_ST(:,dimension)  = obj.A_edgesXj * qF_ST_ij(:,dimension);
            end            
        end
        %%
        function checkIfInDomain(obj)
            if (any(any(obj.Xj)<0) || any(any( ones(obj.N,1)*obj.Omega - obj.Xj <0)))
                error('some points are out of the domain');
            end
        end              
    end
end
