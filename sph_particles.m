%% ToDo:
% - internal energy - non-isothermal
% - boundary conditions (ghost particles)
% - check against c++ code
% - right now: dh/drho = 0
% - regulariseInitialDensity

classdef sph_particles < handle
%% SPH for HVI  - Markus Ganser - TU/e - 2015
% PARTICLES class - defines all function needed to run a sph simulation
    
%% 
    properties        
        N   %amount of particle
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
        %change of density
        drhoj       
        drhoj_diss
        %pressure
        pj             
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
        mj
        %some indicies
        Iin
        Ighost
        Imaterial
        Imaterial_with_boun
        %Material property 
        cj0         
        cj
        beta
        mu
        %scheme
        scheme 
                
        %Time
        t
        dt
        dtfactor %savety factor for timestepping
        tend                
        %Distances rij=xi-xj= rij * xij_h
        rij    % = bar(xij)  
        xij_h  %normalized  = hat(xij)
        nIj
        %kernel
        
        Wij
        Wij_hi
        Wij_hj
        dWij   % to remove later
        dWij_hi
        dWij_hj
        ddWij
        dWij_over_rij
        %smoothing length
        hj 
        h_const  %boolean
        kernel    % M4 | gauss
        %Domain
        Omega
        cell_of_j
        Nc 
        eta2_cutoff        % cut-off-factor
        Rtcell             % size of the cells (first cut-off)
        obsolete_k_pij;
        % boundary condition        
        bc  %struct
        
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

            obj.rho0j = obj_scen.rho0j;            
            obj.rhoj  = obj_scen.rhoj;
            obj.cj0   = obj_scen.c0j;
            obj.cj    = obj_scen.cj;            
            
            obj.beta = obj_scen.beta;
            obj.mu   = obj_scen.mu;
                                    
            obj.mj   = obj_scen.mj;    
                        
            %initial smoothing length
            obj.hj       = ones(obj.N,1)*obj_scen.dx * obj_scen.eta; 
            
            
            %test
%             obj.hj(end-2) = obj.hj(end-2)*1; %2            
%             obj.hj(end-1) = obj.hj(end-1)*1;
        %     obj.hj(end)   = obj.hj(end)*3;%*2.85; %3
        
%              obj.mj(end-2)   = obj_scen.mj(end-2)*1;            ;    
%             obj.mj(end-1)   = obj_scen.mj(end-1)*1;    
%             obj.mj(end)   = obj_scen.mj(end)*0.5;%0.2;    

            
            
            obj.h_const  = obj_scen.h_const;
            obj.kernel   = obj_scen.kernel;
            obj.eta2_cutoff = obj_scen.eta2;
        
            obj.Iin         = obj_scen.Iin;
            obj.Ighost      = [];
            obj.Imaterial   = obj_scen.Imaterial;   
            obj.Imaterial_with_boun= obj.Imaterial;
            obj.bc = obj_scen.bc;     
            
            obj.dtfactor  = obj_scen.dtfactor;
            obj.tend      = obj_scen.tend;           
            obj.vj        = obj_scen.vj;
                                    
            obj.F_int     = zeros(obj.N,obj.dim);
            obj.F_diss    = zeros(obj.N,obj.dim);
            obj.F_diss_art = zeros(obj.N,obj.dim);
            obj.F_ST      = zeros(obj.N,obj.dim);
            obj.nIj       = zeros(obj.N,obj.dim);
            obj.drhoj     = zeros(obj.N,1);
            obj.rhoj_half = zeros(obj.N,1);
            obj.vj_half   = zeros(obj.N,obj.dim);
                  
            obj.pij = [];            
            obj.distances_uptodate = false;
            obj.firststep = true;
            
            obj.cell_of_j      = zeros(obj.N,1);
            obj.obsolete_k_pij = [];

            obj.scheme = obj_scen.scheme;

            obj.t = 0;
            %check
            checkIfInDomain(obj);
            if obj.N > 31000
                warning('quite a lot particles!')
                keyboard
            end
            obj.tmp=obj.mj(end);
        end
        %%
        function start_simulation(obj)            
            disp('#### start of the simulation ####');
            %% time iteration
            icount = 1;
            tic
            ttic = tic;
            obj.IO.initialize()
            while obj.t < obj.tend
                update_dt(obj);
                perform_timestep(obj);
                obj.IO.do(obj);  %plot and save data
                obj.t = obj.t + obj.dt;
                %iterations per seconds
                if toc > 1 %seconds
                    disp(['t=',num2str(obj.t),...
                        's (',num2str(icount),' iter/sec ,Ncon = ',num2str(size(obj.pij,1)),')']);
                    icount = 1;
                    tic    
                    %disp(['tmp: ',num2str(obj.tmp)]);
                end
                icount = icount +1;    
               % keyboard
               % obj.checkDataOnline;
            end
            obj.IO.do(obj);  %plot and save last step
            obj.IO.finalize();
            toc(ttic)
        end
        %%        
        function perform_timestep(obj)
             update_full_step(obj)
             comp_pressure(obj);
             comp_volume(obj);
             
             % mirror step for boundary condition
             if ~isempty(obj.bc)
                mirrorParticles (obj)
             end             
             search_neighbours(obj);
             comp_kernel(obj)             
             comp_forces(obj)
             comp_dRho(obj)                             
             %if ~isempty(obj.damping_area)
             comp_BCdamping(obj)
             %end
             update_half_step(obj)
             update_position(obj)   
             if ~obj.h_const
                 update_h(obj)
             end
%              %% disp something
%              fig=figure(3);
%              NN = 8;
%              plot(obj.Xj(obj.Iin(end-NN:end))',0,'x');
%              
%              hold on;
%              plot(obj.Xj(obj.Iin(end-NN:end))',obj.vj(obj.Iin(end-NN:end)),'x');
%              %diff(obj.Xj(obj.Iin(end-NN:end)))'
%               figpos=fig.Position;
%             figpos(2)=0;
%             fig.Position =figpos;
        end     
        %%
        function update_dt(obj) 
            if obj.dim == 1
                obj.dt = obj.dtfactor * min(obj.hj./(abs(obj.cj + obj.vj)));
            else
                obj.dt = obj.dtfactor * min(obj.hj./(obj.cj + sum(obj.vj(:,1).^2 + obj.vj(:,2).^2,2).^(0.5)));
            end
          %  obj.dt = 1e-5
        end        
        %% neighbour search -  main 
        function search_neighbours(obj)
           cut_off_radius = obj.eta2_cutoff* max(obj.hj);
           if obj.firststep
               obj.Rtcell = cut_off_radius * 1.1; %savety-factor
               disp('create initial cell-structure');               
               initial_search_neighbours(obj);
           elseif (obj.Rtcell < cut_off_radius) %ToDo: make cells smaller
               obj.Rtcell = cut_off_radius * 1.1; %savety-factor
               disp('create new cell-structure');
               initial_search_neighbours(obj);
           else
               update_neighbours(obj);
           end           
        end        
        %%
        function initial_search_neighbours(obj) %upper right and lower left cells are not included
            %% create cell structure
            obj.Nc = floor((obj.Omega(:,2)-obj.Omega(:,1))/obj.Rtcell);
            if obj.dim ==1
                obj.cell_of_j = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                cellshift    = 1;
                Ncell_search = obj.Nc-1;
                if any(obj.cell_of_j ==obj.Nc) ||... %right
                   any(obj.cell_of_j == 1 )           %left
                        warning (' - particels in border cells - algorithm will break - ');
                end                    
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                obj.cell_of_j = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(1));    %convert to 1d counting
                cellshift    = [1,obj.Nc(1)-1,obj.Nc(1),obj.Nc(1)+1]; %upper-left side
                Ncell_search = (obj.Nc(1)*obj.Nc(2) - (obj.Nc(1)+1));
                
                if any(cell_of_j_2d(:,1)==obj.Nc(1)) || ...%in very right cell
                   any(cell_of_j_2d(:,2)==obj.Nc(2)) || ...%in very top cell
                   any(cell_of_j_2d(:,1)==1)         || ...%in very left cell
                   any(cell_of_j_2d(:,2)==1)               %in very bottom cell
                    warning (' - particels in border cells - algorithm will break - ');
                end
                
            end
            %obj.Kcells=sparse(c1d,1:obj.N,1);              %matrix includes a 1-entry if a particle (colum) is in a cell (row)            
            
            if any(obj.cell_of_j > Ncell_search)               
                warning(' - particles not in cell structure! - ')
                keyboard
            end
            
            %% create connectivity list:
            obj.pij = [];            
            sparse_i =[];
            sparse_j =[];
            sparse_s =[];
            sparse_k = 1;
            kp=1;            
            for kcell = 1:Ncell_search %loop over all cells except the upper row                
                j_in_Cell = find(obj.cell_of_j==kcell);
                if size(j_in_Cell,1) > 0   
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
                         
                        %save indice matrices for the connectiviy matrix
                        sparse_i (sparse_k:sparse_k+2*NN-1,:)=...
                                     [i_all(neighbors_of_j+kj);
                                     j_in_Cell(kj)*ones(NN,1)];
                        sparse_j (sparse_k:sparse_k+2*NN-1,:)=...
                                     [(kp:kp+NN-1)';
                                      (kp:kp+NN-1)'];
                        sparse_s (sparse_k:sparse_k+2*NN-1,:)=...
                                     [-1*ones(NN,1);
                                     1*ones(NN,1)];

                        kp=kp+NN;
                        sparse_k=sparse_k+2*NN;
                    end
                end    
            end
            %create connectivity matrix
            obj.AedgesXj =sparse(sparse_i,sparse_j,sparse_s);

            %reset update_neighbours help variables
            obj.obsolete_k_pij = [];
        end    
        %%
        function update_neighbours (obj) 
            %
            function Ashift = lookup_cellshift(obj,jshift)
               if obj.dim==1
                  switch jshift
                      case 1
                          Ashift = 1;
                      case -1
                          Ashift = -1;
                      case inf
                          Ashift = [-1,0,1];
                      otherwise
                            warning('A particel skipped a cell!')
                            keyboard
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
                       case inf      %all surounding cells
                            Ashift = [-1-obj.Nc(1), -obj.Nc(1), 1-obj.Nc(1),...
                                        -1,0,1,...
                                      -1+obj.Nc(1), obj.Nc(1), 1+obj.Nc(1)];
                       otherwise
                            warning('A particel skipped a cell!')
                            keyboard
                   end
               end   
            end
            
            %% create cell structure
            if obj.dim == 1
                cell_of_j_new = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle      
                if any(cell_of_j_new > obj.Nc)  || any(cell_of_j_new < 0)                    
                    warning(' - some particles are not in cell structure anymore! - ')
                    keyboard
                end
           
            else
                cell_of_j_2d = floor(ones(obj.N,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj - ones(obj.N,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                cell_of_j_new = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(1));    %convert to 1d counting
                if any(cell_of_j_2d(1) > obj.Nc(1) || cell_of_j_2d(2) > obj.Nc(2))...
                              || any(any(cell_of_j_2d < 0))                    
                    warning(' - some particles are not in cell structure anymore! - ')
                    keyboard
                end
            end
                                                        
            
            %%
            np_change = size(cell_of_j_new,1)-size(obj.cell_of_j,1); %
 
            % adjust vector sizes
            if np_change > 0 %new particles
                obj.cell_of_j = [obj.cell_of_j; inf*ones(np_change,1)];
                %enlarge matrix in order to have as much rows as particles
                obj.AedgesXj(obj.N,1) = 0;   
            elseif np_change <0 %particles to remove
                cell_of_j_new = [cell_of_j_new; inf*ones(-np_change,1)];
            end
            cell_shift_of_j = cell_of_j_new - obj.cell_of_j;
            
            if any(isnan(cell_shift_of_j))
                keyboard
            end
            
            j_changes_cell = find(cell_shift_of_j); %alternative: a-b ~=0
                              
            %% loop over all particels which have changed the cell
            for j = j_changes_cell'  
                
                if cell_shift_of_j(j) == -inf % new particle
                    jghost = 1; %(1: new, 0-moved, -1: removed)
                    jreal = false;
                    cell_shift_of_j(j) = inf;
                elseif  cell_shift_of_j(j) == inf %obsolete particle
                    jghost = -1;
                    jreal = false;
                    cell_shift_of_j(j) = inf;
                else
                    if any(j==obj.Ighost);
                        jghost = 0;
                        cell_shift_of_j(j) = inf;
                        jreal = false;
                    else
                        jreal = true;
                        jghost = nan;
                    end
                end
                
                cellshift = lookup_cellshift(obj,cell_shift_of_j(j)); %cells ahead in moving direction (for regular particles)
                
                %% find new and obsolete particles:
                
                %find obsolete points 
                if (jghost<=0 || jreal)
                    if jreal %normal particles
                        % (all particels in cells which are no
                        %neighbours any more - right now, one behind ...|X|j|o|o|...
                        cells_to_look_at = obj.cell_of_j(j) - cellshift;
                    else %remove all old connectivities, cellshift are all sourinding cells
                        cells_to_look_at = obj.cell_of_j(j) + cellshift; %
                    end
                    i_in_obs_neighb_cell = (double(ismember(obj.cell_of_j,cells_to_look_at))); 

                    %find the belonging connectivities                    
                    i_in_obs_neighb_cell(j,1) = 5;  % just a number,to obtain an unsymmetric connection - abs(+-5 -+1)=4
                    temp    = (i_in_obs_neighb_cell'*obj.AedgesXj);                    
                    pij_obs = (abs(temp)==4);
                    new_obsolete_k_pij = find(pij_obs');
                else
                    new_obsolete_k_pij = [];
                end
                % move j to new cell (each particel at one time for
                % consistency)     % ...|o|j|o|o|... ->  ...|o|o|j|o|...                
                obj.cell_of_j(j) = cell_of_j_new(j);               
                
                if (jghost>=0 || jreal)
                    %find new points in cells ahead ...|o|o|j|X|...
                    cells_to_look_at = obj.cell_of_j(j) + cellshift;
                    i_in_new_neighb_cell = find(ismember(obj.cell_of_j,cells_to_look_at));
                    if ~jreal
                        % in case of a ghost particle, remove the
                        % selfconnection (not necessary for real particle
                        % since I look only in neigbouring cells)
                        i_in_new_neighb_cell=...
                            i_in_new_neighb_cell(i_in_new_neighb_cell~=j);
                    end
                    nCon_new = size(i_in_new_neighb_cell,1);     
                else
                    nCon_new = 0;
                end
                                
                %% apply changes %%
                
                % clear all obsolete connectivities of node j in the
                % adjacency matrix nodes <-> edges
                obj.AedgesXj(:,new_obsolete_k_pij) = 0;
                %set a self connection for obsoltete particles
                obj.pij(new_obsolete_k_pij,:) = ...
                             ones(size(new_obsolete_k_pij,1) ,1) *[1,1]; 
                obj.obsolete_k_pij = [obj.obsolete_k_pij;...
                                         new_obsolete_k_pij];
                % http://de.mathworks.com/matlabcentral/newsreader/view_thread/283923

                % insert new connectivies               
                if nCon_new > 0
                    pij_new = [j*ones(nCon_new,1),i_in_new_neighb_cell]; 
                    nCon_obsolete = length(obj.obsolete_k_pij);
                    if nCon_new <= nCon_obsolete
                         k_pij_new =  obj.obsolete_k_pij(1:nCon_new);
                         obj.pij(k_pij_new,:) = pij_new;   
                         obj.obsolete_k_pij = obj.obsolete_k_pij(nCon_new+1:end);                            
                    else
                         %write on present connectivities
                         obj.pij(obj.obsolete_k_pij,:)= pij_new (1:nCon_obsolete,:);
                         %and create new ones
                         n_pij_new = nCon_new - nCon_obsolete;
                         obj.pij   = [obj.pij;
                                     pij_new(nCon_obsolete+1:nCon_new,:)];
                         n_Apij    = size(obj.pij,1);
                         k_pij_new = [obj.obsolete_k_pij;...
                                     ((n_Apij-n_pij_new+1):n_Apij)'];
                                 
                         % enlarge adjacency matrix 
                         obj.AedgesXj(obj.N,n_Apij) = 0;   
                         
                         % remove obsolete flags
                         obj.obsolete_k_pij = []; 
                    end

                    %% update adjacency matrix                    
                    % write new connectivities in adjacency matrix
                    obj.AedgesXj (i_in_new_neighb_cell,k_pij_new) = -1*speye(nCon_new); %ingoing  
                    obj.AedgesXj (j , k_pij_new) = 1; %outgoing                    
                end
                
            end 
            
            % remove the obsolete ghost particles (which are =inf) in cell_of_j;
            if np_change < 0               
               obj.cell_of_j = obj.cell_of_j(1:end-(-np_change));
               %and shrink connectivity matrix          
               obj.AedgesXj(end-(-np_change)+1:end,:)=[];  
            end
        end                        
        %%             
        function comp_distances(obj)
            obj.xij_h = (obj.Xj(obj.pij(:,1),:)-obj.Xj(obj.pij(:,2),:));
            obj.rij   = sum(obj.xij_h.^2,2).^0.5;
            obj.xij_h = obj.xij_h./(obj.rij*ones(1,obj.dim));  
            
            % flag particle which are within the cutoff radius           ^       
            temp = (obj.rij < obj.eta2_cutoff *...
                      max(obj.hj(obj.pij(:,1)),obj.hj(obj.pij(:,2))))...
                   .*(obj.rij > 0);
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
        %% compute kernel values incl. derivatives
        function comp_kernel(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
             if obj.h_const && strcmp(obj.kernel,'Gauss')% if h is constant
                comp_kernel_gauss_hconst(obj);
                obj.Wij_hi = obj.Wij;
                obj.Wij_hj = obj.Wij;
                obj.dWij_hi = obj.dWij;
                obj.dWij_hj = obj.dWij;
            else
                if strcmp(obj.kernel,'Gauss')
                    sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
                    fw  = @(r) exp(-r.^2);
                    fdw = @(r)  -2*r.*exp(-r.^2);
                elseif strcmp(obj.kernel,'M4')
                    sigma = [2/3; 10/(7*pi); 1/pi];
                    fw = @(r) ...
                        1/4*(2-r).^3   ...
                         - (1-r).^3 .* (r<1);            
                    fdw = @(r) ... 
                         -3/4*(2-r).^2  ...
                        + 3*(1-r).^2 .* (r<1);  
                elseif strcmp(obj.kernel,'Wendland')
                    sigma = [3/4; 7/(4*pi); 21/(16*pi)];
                    fw = @(r)  ...
                       (1-r/2).^4 .* (1+2*r);
                    fdw =@(r)... 
                      -5*r.*(1-r/2).^3;    
                else
                    error([obj.kernel,' as kernel is not implemented yet']);
                end
                
                r = [(obj.rij(obj.active_k_pij)) ./ obj.hj(obj.Ii);
                     (obj.rij(obj.active_k_pij)) ./ obj.hj(obj.Ij)];
 
                Np = sum(obj.active_k_pij);
                I = r < 2;
                w = zeros(2*Np,1);            
                w(I) = sigma(obj.dim)*fw(r(I));

                obj.Wij_hi = w(1:Np)./ (obj.hj(obj.Ii).^obj.dim);
                obj.Wij_hj = w(Np+1:end)./ (obj.hj(obj.Ij).^obj.dim);
                obj.Wij    = 0.5*(obj.Wij_hi + obj.Wij_hj);

                dw = zeros(2*Np,1);
                dw(I,:)   = sigma(obj.dim)*fdw(r(I,:));

                obj.dWij_hi = 1./obj.hj(obj.Ii).^(obj.dim+1) .* dw(1:Np);
                obj.dWij_hj = 1./obj.hj(obj.Ij).^(obj.dim+1) .* dw(Np+1:end);
                obj.dWij    = 0.5*(obj.dWij_hi+obj.dWij_hj);
                
                if obj.beta ~= 0 % necessary only for surface tension
                    error('ddW not implented for not constant h yet');
%                 obj.dWij_over_rij = obj.dWij ./ obj.rij(obj.active_k_pij); 
% 
%                 obj.ddWij = 1/obj.hj^2 * sigma(obj.dim) .* ...
%                     (r<2)  .*(...
%                     +6/4*(2-r)  ...
%                     - 6*(1-r) .* (r<1)...
%                     ); 
                end   
             end
        end         
        %% Gauss kernel
        function comp_kernel_gauss_hconst(obj)
            sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
            obj.Wij   = 1/(obj.hj(1)^obj.dim)* sigma(obj.dim) *...
                        exp(-(obj.rij(obj.active_k_pij)/obj.hj(1)).^2);
            obj.dWij  = obj.Wij.* -2.*obj.rij(obj.active_k_pij)/obj.hj(1)^2;
            obj.dWij_over_rij = obj.Wij.* -2./obj.hj(1)^2;
            obj.ddWij = -2/obj.hj(1)^2*(obj.dWij.* obj.rij(obj.active_k_pij) + obj.Wij );
        end       
        %%    
        function update_full_step(obj) 
            if ~obj.firststep
                obj.rhoj(obj.Iin) = obj.rhoj_half(obj.Iin)+ 0.5*obj.dt * obj.drhoj(obj.Iin);
                obj.vj(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + 0.5*obj.dt *...
                    (obj.F_int(obj.Iin,:)+obj.F_ST(obj.Iin,:)+obj.F_diss(obj.Iin,:))...
                    ./(obj.mj(obj.Iin)*ones(1,obj.dim));
            end
        end
        %%
        function update_half_step(obj) 
            if obj.firststep
                obj.rhoj_half(obj.Iin) = obj.rhoj(obj.Iin)  + 0.5*obj.dt*obj.drhoj(obj.Iin);
                obj.vj_half(obj.Iin,:)  = obj.vj(obj.Iin,:) + 0.5*obj.dt * ...
                    (obj.F_total(obj.Iin,:))./(obj.mj(obj.Iin)*ones(1,obj.dim));
                obj.firststep=false;
            else
                obj.rhoj_half(obj.Iin) = obj.rhoj_half(obj.Iin)+ obj.dt * obj.drhoj(obj.Iin);
                obj.vj_half(obj.Iin,:)  = obj.vj_half(obj.Iin,:) + obj.dt *...
                    (obj.F_total(obj.Iin,:))./(obj.mj(obj.Iin)*ones(1,obj.dim));
                %+ones(size(Iin))*g_ext);
            end
        end
        %%
        function update_position(obj) 
             obj.Xj(obj.Iin,:)      = obj.Xj(obj.Iin,:) + obj.dt * obj.vj_half(obj.Iin,:);
             obj.distances_uptodate = false;
        end
        %%
        function update_h(obj)
           obj.hj(obj.Iin) = obj.hj(obj.Iin)...
               - obj.dt * (obj.hj(obj.Iin)./obj.rhoj(obj.Iin)) .* obj.drhoj(obj.Iin);
           %disp(['min: ',num2str(min(obj.hj)),'; max: ',num2str(max(obj.hj))]);
        end
        %%
        function comp_volume(obj)
            %only necessary for v-scheme
            obj.Vj = obj.mj(:)./obj.rhoj(:);
        end        
        %%
        function comp_pressure(obj) %switch
            comp_pressure_isothermal(obj)
        end
        %%
        function comp_pressure_isothermal(obj)
            obj.pj(obj.Iin,:) = obj.cj0(obj.Iin) .* obj.cj0(obj.Iin) ...
                                 .*(obj.rhoj(obj.Iin) - obj.rho0j(obj.Iin));
        end
        %%
        function comp_dRho(obj)
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            if strcmp(obj.scheme,'m')            
                comp_dRho_m_scheme(obj);
            elseif strcmp(obj.scheme,'v')
                comp_dRho_v_scheme(obj);
            else
                error('scheme not supported');
            end
            
            % compute dissipative mass flux
            comp_dRho_diss(obj)
            obj.drhoj = obj.drhoj + obj.drhoj_diss;          
        end
        %%
        function comp_dRho_m_scheme(obj)
            obj.drhoj  = 0*obj.drhoj; 
            qrho_ij = zeros(size(obj.pij,1),1);
            qrho_ji = zeros(size(obj.pij,1),1);
            qrho_ij(obj.active_k_pij,:) =  ... % i-> j %hi
                obj.mj(obj.Ij) .*... 
                sum(... %scalar product v*n 
                (obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)) .*...
                (obj.dWij_hi*ones(1,obj.dim)).*obj.xij_h((obj.active_k_pij),:)...
                ,2);
            qrho_ji(obj.active_k_pij,:) =  ... % j -> i %hj
                obj.mj(obj.Ii).*...
                 sum(... %scalar product v*n 
                (obj.vj(obj.Ij,:)-obj.vj(obj.Ii,:)) .*...
                (obj.dWij_hj*ones(1,obj.dim)).*-obj.xij_h((obj.active_k_pij),:)...
                ,2);
                        
            %add up all the corresponding density flux in each node
            obj.drhoj = ((obj.AedgesXj>0) * qrho_ij +...
                        ( obj.AedgesXj<0) * qrho_ji);  %ToDo: make more computational effective
        end
        %%
        function comp_dRho_v_scheme(obj)
            %todo: differnt dWij (for hj and hi)
            obj.drhoj  = 0*obj.drhoj; 
            qrho_ij = zeros(size(obj.pij,1),1);
            qrho_ji = zeros(size(obj.pij,1),1);
            qrho_ij(obj.active_k_pij,:) =  ... % i-> j %hi
                obj.Vj(obj.Ij) .*... 
                sum(... %scalar product
                (obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)) .*...
                (obj.dWij_hi*ones(1,obj.dim)).*obj.xij_h((obj.active_k_pij),:)...
                ,2);
            qrho_ji(obj.active_k_pij,:) =  ... % j -> i %hj
                obj.Vj(obj.Ii).*...
                 sum(... %scalar product 
                (obj.vj(obj.Ij,:)-obj.vj(obj.Ii,:)) .*...
                (obj.dWij_hj*ones(1,obj.dim)).*-obj.xij_h((obj.active_k_pij),:)...
                ,2);
                        
            %add up all the corresponding density flux in each node
            obj.drhoj = ((obj.AedgesXj>0) * qrho_ij +...
                         (obj.AedgesXj<0) * qrho_ji)...
                        .*obj.rhoj/1;
        end
        %%
        function comp_dRho_diss(obj) % Dissipative mass flux (Zizis2014)
            alpha_mass = 0.3;  %super sensitive!
            obj.drhoj_diss  = 0*obj.drhoj_diss; 
            qrho_ij_diss = zeros(size(obj.pij,1),1);
            qrho_ji_diss = zeros(size(obj.pij,1),1);
            
            rho_ij_bar = (obj.rhoj(obj.Ii) + obj.rhoj(obj.Ij))/2;
            %approach 1 (remove vijxijh<=0 in qrho_ij_diss!)
            cj_bar     = (obj.cj(obj.Ii)   + obj.cj(obj.Ij)  )/2;
            v_sig      = cj_bar;
            %approach 2 (switch on vijxijh<=0 in qrho_ij_diss!)
%             beta_mass  = 1;
%             vijxijh = sum(...
%                    ((obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:))).*... (vi-vj)*hat(xij)
%                      obj.xij_h(obj.active_k_pij,:)...
%                      ,2);
%              v_sig1 = (0.5 *(obj.cj(obj.Ii)+obj.cj(obj.Ij))...
%                    - 0.5*beta_mass*vijxijh); 
% 
%             
            qrho_ij_diss(obj.active_k_pij,:) = alpha_mass* ... % i-> j %hi
                obj.mj(obj.Ij).*...
                (obj.pj(obj.Ii)-obj.pj(obj.Ij))./...
                (rho_ij_bar .* v_sig)...
                .*obj.dWij;%.*(vijxijh<=0);
           
            qrho_ji_diss(obj.active_k_pij,:) = alpha_mass* ... % i-> j %hi
                obj.mj(obj.Ii).*...
                (obj.pj(obj.Ij)-obj.pj(obj.Ii))./...
                (rho_ij_bar .* v_sig)...
                .*obj.dWij;%.*(vijxijh<=0);            
            % spread flux into nodes
            obj.drhoj_diss = ((obj.AedgesXj>0) * qrho_ij_diss +...
                              (obj.AedgesXj<0) * qrho_ji_diss);
                      
%             if any(obj.drhoj_diss > 1e3)
%                 keyboard
%             end
        end
        %%
        function comp_forces(obj)
            %internal forces
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            if strcmp(obj.scheme,'m')            
                comp_Fint_m_scheme(obj);
            elseif strcmp(obj.scheme,'v')
                comp_Fint_v_scheme(obj);
            else
                error('scheme not supported');
            end
            
            %artifisial dissipation
            comp_Fdiss_art(obj)
            
            %dissipation
            if obj.mu ~=0
                comp_Fdiss(obj);
            end
            % surface tension
            if obj.beta ~= 0
                comp_F_ST(obj);
            end
            %add up everything
            obj.F_total(obj.Iin,:) = obj.F_int(obj.Iin,:)...
                                 + obj.F_diss(obj.Iin,:)...
                                 + obj.F_diss_art(obj.Iin,:)...
                                 + obj.F_ST(obj.Iin,:);
%             if any(sum(obj.F_total)>1e-10)
%                 keyboard
%             end
        end
        %%
        function comp_Fint_m_scheme(obj) %internal forces  
            %compute the flux
            qF_int_ij = zeros(size(obj.pij,1),obj.dim);
            p_rho_omega_j = obj.pj./(1*obj.rhoj.*obj.rhoj);
            
            qF_int_ij(obj.active_k_pij,:)=...
                -((obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... % hi
                (p_rho_omega_j(obj.Ii).*obj.dWij_hi)...  %<- change here for variable h
                *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:)...
                +...
                (obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... %hj
                (p_rho_omega_j(obj.Ij).*obj.dWij_hj)...  %<- change here
                *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:));

            %spread fluxes to the nodes
            obj.F_int = obj.AedgesXj * qF_int_ij;  
        end
        %%
        function comp_Fint_v_scheme(obj) %internal forces  
            %compute the flux
            qF_int_ij = zeros(size(obj.pij,1),obj.dim);
            p_omega_j = obj.pj./(1);
            
            qF_int_ij(obj.active_k_pij,:)=...
                -((obj.Vj(obj.Ii).*obj.Vj(obj.Ij) .*... % hi
                (p_omega_j(obj.Ii).*obj.dWij_hi)...  %<- change here
                *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:)...
                +...
                (obj.Vj(obj.Ii).*obj.Vj(obj.Ij) .*... %hj
                (p_omega_j(obj.Ij).*obj.dWij_hj)...  %<- change here
                *ones(1,obj.dim)).*obj.xij_h(obj.active_k_pij,:));

            %spread fluxes to the nodes
            obj.F_int = obj.AedgesXj * qF_int_ij;
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
            obj.F_diss = obj.AedgesXj * qF_diss_ij;

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
                 ((obj.mj(obj.Ii) .* obj.mj(obj.Ij).* ...
                 alpha_diss .*v_sig .* vijxijh) ./...
                 (0.5 *(obj.rhoj(obj.Ii) + obj.rhoj(obj.Ij)))... bar(rho)
                 .*obj.dWij)...   
                 *ones(1,obj.dim)...
                 .*obj.xij_h(obj.active_k_pij,:);
             
            obj.F_diss_art = obj.AedgesXj * qF_diss_art_ij;
        end        
        %%
        function comp_F_ST(obj)  %Violeau423        
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
            obj.F_ST  = obj.AedgesXj * qF_ST_ij;         
        end
        %%
        function comp_BCdamping(obj) %todo
           for boun = obj.bc
               if ~isempty(boun.damping_area)
%                    xb = obj.bc.damping_area;
%                    L  = xb(2)-xb(1);           
%                    kb = logical((obj.Xj(obj.Iin)>xb(1)) .* (obj.Xj(obj.Iin)<xb(2))); %particels in boundary layer
%                    x  = obj.Xj(kb,:); 
%                    Nkb = sum(kb);   
%                    %hyperbolic
%                    sigma0 = max(obj.cj) / L;
%                    sigma  = sigma0 *(x - xb(1)) ./ (xb(2)-x +0.5*max(obj.hj));
%         %             
%         %            m      = 10;
%         %            sigma0 = 10;
%         %            sigma=sigma0*((obj.Xj(kb,:) - xb(1))./(xb(2)-xb(1))).^m;
% 
%                     %switching functiosn
%                    f0 = ones(Nkb,1);             % do nothing
%                    f1 = - (x  - xb(2)) / L;         % linear
%                    f2 = (L^2 - (x -xb(1)).^2) / L^2;  % quadratic
% 
%                       %% damp if
%                       %never
%                    damp_never = zeros(Nkb,1);
%                       %always           
%                    damp_always = ones(Nkb,1);
%                       %negative pressure
%                    damp_pressure = obj.pj(kb,:)<0; 
%                       %only if v<0
%                    damp_velocity = obj.vj(kb,:)<0; 
% 
%                    %% setup
%                    damp = damp_never;
%                    f_S = f1;  %density-change
%                    f_Q = f1;  %forces
%                    sigmaS = 1*sigma;
%                    sigmaQ = 0.02*sigma;
% 
% 
%                     %density:
%                    S = - sigmaS.*(obj.rhoj(kb,:) - obj.rho0j(kb,:));
%                    obj.drhoj(kb,:) = (obj.drhoj(kb,:).*f_S + S.*damp);         
%                     %forces:
%                    Q = - sigmaQ.*(obj.vj(kb,:) - 0);
%                    obj.F_total(kb,:) = (obj.F_total(kb,:).*f_Q + Q.*damp ) ;
% 
% 
%                    %% stop particels
%         %            kb_stop = logical ((obj.Xj>=xb(2)));
%         %            obj.vj(kb_stop) = 0;
               end
           end
        end       
        %%                       
        function mirrorParticles (obj) %todo
   
    %todo: use cell-structure:
%                     cells_to_look_at = obj.cell_of_j(obj.mirrorParticlesj)
%                     j_to_lock_at =...
%                         (double(ismember(obj.cell_of_j,cells_to_look_at))); 
%             if obj.dim == 1
%                 Ashift = [-1,0,1];
%             else
%                 Ashift = [-1-obj.Nc(1), -obj.Nc(1), 1-obj.Nc(1),...
%                                -1,0,1,...
%                           -1+obj.Nc(1), obj.Nc(1), 1+obj.Nc(1)];
%             end
            
            function d = distanceToParticles(x,xborder,outer_normal)
                borderline = mean(xborder); 
                if size(x,2) == 1
                    d = abs( x - borderline);
                else
                    if outer_normal(1)~=0 && outer_normal(2) ==0
                        k=1; %vertical
                    elseif outer_normal(1)==0 && outer_normal(2) ~=0
                        k=2; %horizontal
                    else
                        error('other normal-directions are not implemented yet');
                    end
                    d = abs( x(:,k) - borderline(k));
                end
            end
            
 
            function d = point_to_line(pt, v1, v2)
                  Np = size(pt,1);
                  a = [v1 - v2,0];
                  b = [pt - ones(Np,1)*v2,zeros(Np,1)];
                  d = sum(cross(ones(Np,1)*a,b).^2,2).^0.5 ./ norm(a);  
            end
            
           % keyboard
            %reset to initial
            obj.Ighost = [];
            obj.Xj   = obj.Xj(obj.Iin,:);
            obj.vj   = obj.vj(obj.Iin,:);
            obj.mj   = obj.mj(obj.Iin);
            obj.Vj   = obj.Vj(obj.Iin);
            obj.rhoj = obj.rhoj(obj.Iin);
            obj.hj   = obj.hj(obj.Iin);
            obj.cj   = obj.cj(obj.Iin);
            obj.pj   = obj.pj(obj.Iin);
            obj.N    = size(obj.Xj,1);
            obj.Imaterial_with_boun = obj.Imaterial;  
            
            %add boundaries:            
            kk=0;
            for boun = obj.bc
                if ~isempty(boun.mirrorParticlesj) && isempty(boun.p1) %nr                                                        
                    kb = boun.mirrorParticlesj;
                    ki = find(kb==0);  %ki = obj.Iin;
                    outer_normal = boun.outer_normal;
                    
                    if obj.dim == 1
                        v_factor = -1;
                    else
                        v_factor = [-1,-1];
                    end
                    p_factor = -1   ;
                    
                    j_to_lock_at = obj.Iin(ki); %all particles
                    % can probably be improved when taking the cell structure into
                    % account
                    d = distanceToParticles(obj.Xj(j_to_lock_at,:),...
                                            obj.Xj(obj.Iin(kb),:),...
                                            outer_normal);

                    Imirror_local = d <  obj.eta2_cutoff * max(obj.hj(obj.Iin));% * 12;
                    Imirror       = obj.Iin(ki(Imirror_local));
                    
                     % to examine:
                    if obj.vj(obj.Iin(end),1)>0  
                        mV_factor_min = 0.3;
                        mV_factor_max = 1.1;
                        dmax = max(d(Imirror_local));
                        drel = (dmax - d(Imirror_local))/dmax;
                        mV_factor = (1-drel)*mV_factor_min +...
                                     drel *mV_factor_max;
                        %mV_factor = ones(size(Imirror,1),1);
                    else
                        mV_factor_min = 2;
                        mV_factor_max = 0.2;
                        dmax = max(d(Imirror_local));
                        drel = (dmax - d(Imirror_local))/dmax;
                        mV_factor = (1-drel)*mV_factor_min +...
                                     drel *mV_factor_max;
                        %mV_factor = ones(size(Imirror,1),1);
                    end
                elseif isempty(boun.mirrorParticlesj) && ~isempty(boun.p1)
                    outer_normal = boun.outer_normal;
                    if obj.dim ==1
                        v_factor = -1;
                    else
                       if outer_normal(1)~=0 && outer_normal(2) ==0
                            v_factor = [-1,1];
                       elseif outer_normal(1)==0 && outer_normal(2) ~=0
                            v_factor = [1,-1];
                       else
                          error('other normal-directions are not implemented yet');
                       end
                    end
                    
                    mV_factor  = 1;
                    p_factor = 1;
                    j_to_lock_at = obj.Iin; %all particles
                    pXj = obj.Xj(j_to_lock_at,:);
                    if obj.dim == 1
                         d = abs(pXj-boun.p1);
                    else
                         d = point_to_line(pXj,boun.p1,boun.p2);
                    end
                    Imirror_local = d <  obj.eta2_cutoff * max(obj.hj(obj.Iin));% * 12;
                    Imirror       = obj.Iin(Imirror_local);
                else
                    error('no such boundary condition possible');
                end
                NN = size(Imirror,1);                   
                obj.Ighost(kk+(1:NN),:) = Imirror;
                
                %copy data
                if ~isempty(Imirror)
                    obj.Xj (obj.N+(1:NN),:) = ...
                             obj.Xj(Imirror,:) + 2*d(Imirror_local)*outer_normal; 

                    obj.vj (obj.N+(1:NN),:) = ...
                              obj.vj(Imirror,:).*(ones(NN,1)*v_factor);  %probably only valid in 1D!
                    obj.pj (obj.N+(1:NN),:) = ...
                              obj.pj(Imirror).*p_factor; 

                    %to examine:
                    mirror_mass=obj.mj(Imirror);
                    %mirror_mass=obj.tmp;
                    obj.mj (obj.N+(1:NN),:) = ...
                              mirror_mass.*mV_factor;

                    obj.Vj (obj.N+(1:NN),:) = ...
                              obj.Vj(Imirror).*mV_factor;

                    obj.rhoj (obj.N+(1:NN),:) = ...
                                obj.rhoj(Imirror);
                    obj.hj (obj.N+(1:NN),:) = ...
                              obj.hj(Imirror);
                    obj.cj (obj.N+(1:NN),:) = ...
                              obj.cj(Imirror); 

                    obj.N = size(obj.Xj,1);

                    obj.Imaterial_with_boun=[obj.Imaterial_with_boun;
                                             obj.N-NN+1 , obj.N];
                end
                kk=kk+NN;
            end
        end            
        %%
        function checkIfInDomain(obj)
            if (any(any(obj.Xj-  ones(obj.N,1)* obj.Omega(:,1)'<0)) || any(any( ones(obj.N,1)*obj.Omega(:,2)' - obj.Xj <0)))
                error('some points are out of the domain');
            end
        end  
        %%
        function checkDataOnline(obj)                       
            disp ('check data');
            %%
            checkIfInDomain(obj)
            %%
            %check if AedgesXj has equally many in and outgoing
            %connectivities
            if any(find(sum(obj.AedgesXj,1)))
                keyboard
            end
            %%
        end
    end
end
