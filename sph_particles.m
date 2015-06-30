%% ToDo:
% - check against c++ code
% - regulariseInitialDensity
% compGhost is necessary - otherwise noflow-bc breaks down - why?

%--
% - axisymetrics
% - shear stresses
% - damage model
% - thermal diffusion? heat diffuses in solids quite fast... what is that
% timescale?


classdef sph_particles < handle
%% SPH for HVI  - Markus Ganser - TU/e - 2015
% PARTICLES class - defines all function needed to run a sph simulation
    
%% 
    properties        
        N   %amount of particle
        dim
        %forces
        Fj_tot
        Fj_int
        Fj_diss
        Fj_ST
        Fj_phy_diss
        %density
        rho0j
        rhoj               
        rhoj_half
        rhoj_real  %for axissymetric use
        %change of density
        drhoj_tot
        drhoj_int       
        drhoj_diss
        %pressure
        pj             
        %velocity
        vj
        vj_half
        %position
        Xj
        %internal energy
        ej
        ej_half
        % change of energy
        dej
        dej_diss
        %
        c1j_half
        c2j_half
        c3j_half
        c4j_half
        c4j
        dc4j
        %
        %connectivity ([indice of startpoint, indice of endpoint]^ammount of connections
        pij  
        %all active connection of pij in two seperate arrays 
        active_k_pij % all edges inside the cutoffradius
        Ii    %indice of point i of the ij-edge
        Ij    %indice of point j of the ij-edge
        %Incidence matrix: row: particles, column: edges (with respect to A_pij)
        AedgesXj;  
        %volume
        Vj
        %mass
        mj
        %some indicies
        Iin
        Ighost
        Icomp %what particles to compute
        Imaterial
        Imaterial_with_boun
        %Material property 
        c0j         
        cj
        MG_Gammaj   %gruneisen parameter
        MG_Sj       % relation between shock and particle velocity
        
        beta
        mu
        art_diss_para
        
        g_ext    %gravity

        %scheme
        scheme 
        % equation of state
        EOS
        %Time
        t
        dt
        dt_fix
        dtfactor %savety factor for timestepping
        tend  
        tpause
        %Distances rij=xi-xj= rij * xij_h
        rij    % = bar(xij)  
        xij_h  %normalized  = hat(xij)
        nIj
        %kernel
        
        fw   %kernel function handle
        fdw  %kernel function handle
        Wij
        Wij_hi
        Wij_hj
        dWij   % to remove later
        dWij_hi
        dWij_hj
        ddWij
        dWij_over_rij
        %smoothing length
        eta
        hj 
        Omegaj
        compOmegaj %bool
        h_const  %boolean
        kernel    % M4 | gauss
        kernel_cutoff
        %Domain
        Omega
        cell_of_j
        Nc 
        cell_cutoff        % cut-off-factor
        Rtcell             % size of the cells (first cut-off)
        obsolete_k_pij;
        % boundary condition        
        bc  %struct
        %        
        % experimantal settings (struct)
        exp_settings
        Xj_ghost
        pij_ghost
        %some flags
        isothermal 
        distances_uptodate
        firststep
          % compute properties of the ghost varialbes
            %(although they are getting deleted )
            % needed to compute the extraction of momentum/ernergy
        compGhost
        % IO class
        IO    
        
        %axis-symmetric variables:
        flag_axi
        f1 %function
        dxf1 %function
        f1j  %computed values
        dxf1j
        %
        tmp
        IOdata %(pj(Ighost),...) struct
        IOdata_h % halfstep data
    end
 %%   
    methods
       
        %% constructor
        function obj = sph_particles(obj_scen)
            if nargin==0
                disp('##########################');
                disp('#### start 1d_riemann ####');
                disp('##########################');
                scen_2d_impact_axi();
                return
            end
            %initialize IO and read data if necessary            
            obj.IO = sph_IO(obj_scen); 
           
            obj.Omega = obj_scen.Omega;
            obj.Xj    = obj_scen.Xj;
            obj.N     = size(obj.Xj,1);
            obj.dim   = size(obj.Xj,2);

            obj.rho0j = obj_scen.rho0j;            
            obj.rhoj  = obj_scen.rhoj;
            obj.c0j   = obj_scen.c0j;
            obj.cj    = obj_scen.cj;            
            
            obj.beta  = obj_scen.beta;
            obj.mu    = obj_scen.mu;
            obj.art_diss_para = obj_scen.art_diss_para;
            obj.g_ext = obj_scen.g_ext;                      
            obj.mj    = obj_scen.Vj.*obj_scen.rhoj;    
            obj.Vj    = obj_scen.Vj;            
            %initial smoothing length
            obj.eta      = obj_scen.eta;
            obj.hj       = obj_scen.Vj.^(1/obj.dim) * obj.eta; 
            obj.Omegaj   = 1; % let it scalar if not used
            obj.compOmegaj = obj_scen.compOmegaj;

            obj.exp_settings = obj_scen.exp_settings;                                 
            
            obj.h_const       = obj_scen.h_const;
            obj.kernel        = obj_scen.kernel;
            obj.kernel_cutoff = obj_scen.kernel_cutoff;
        
            obj.Iin         = obj_scen.Iin;
            obj.Icomp       = obj.Iin;
            obj.Ighost      = [];
            obj.Imaterial   = obj_scen.Imaterial;   
            obj.Imaterial_with_boun= obj.Imaterial;
            
            obj.bc        = obj_scen.bc;     
           
            obj.dt        = obj_scen.dt;
            obj.dt_fix    = obj_scen.dt;
            obj.dtfactor  = obj_scen.dtfactor;
            obj.tend      = obj_scen.tend;           
            obj.vj        = obj_scen.vj;
            obj.tpause    = obj_scen.tpause;           
            obj.ej        = obj_scen.ej;
            obj.MG_Gammaj = obj_scen.MG_Gammaj;
            obj.MG_Sj     = obj_scen.MG_Sj;

            obj.scheme = obj_scen.scheme;
            
            if strcmp(obj.scheme,'a')
                obj.flag_axi = true;
            else
                obj.flag_axi = false;
            end
            obj.EOS    = obj_scen.EOS;
            
            %compute the energy only when needed:
            obj.isothermal=false;
          %  obj.isothermal = strcmp(obj.EOS,'ISO') || strcmp(obj.EOS,'Water');
            
         
            %define kernel:
            if strcmp(obj.kernel,'Gauss')
                sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
                obj.fw  = @(q) sigma(obj.dim) * exp(-q.^2);
                obj.fdw = @(q) sigma(obj.dim) * -2*q.*exp(-q.^2);
            elseif strcmp(obj.kernel,'M3')
                sigma = [2/3; 10/(7*pi); 1/pi];
                obj.fw = @(q)  sigma(obj.dim)* (...
                    1/4*(2-q).^3   ...
                      -( (1-q).^3 .* (q<1)));            
                obj.fdw = @(q)  sigma(obj.dim)* (... 
                     -3/4*(2-q).^2  ...
                    + 3*(1-q).^2 .* (q<1));  
            elseif strcmp(obj.kernel,'M4')
                sigma = [1/24, 96/(1199*pi), 1/(20*pi)];
                obj.fw = @(q)  sigma(obj.dim)* (...
                    (5/2-q).^4 ...
                    -5*(3/2-q).^4 .* (q<1.5) ...
                    +10*(0.5-q).^4 .*(q<0.5));
                obj.fdw = @(q)  sigma(obj.dim)* -4* (... 
                    (5/2-q).^3 ...
                   -5 *(3/2-q).^3 .* (q<1.5) ...
                   +10*(0.5-q).^3 .* (q<0.5));      
               
               %axisymmetric 
                obj.f1 =@(xi) ...
                    (0<=xi).*(xi<1).*...
                        (7/15*xi.^-1 + 2/3*xi - 1/6*xi.^3 + 1/20*xi.^4).^-1 +...
                    (1<=xi).*(xi<2).*...
                        (8/15*xi.^-1 - 1/3 + 4/3*xi - 2/3*xi.^2 + 1/6*xi.^3 - 1/60*xi.^4).^-1 +...
                    (2<=xi).*1;            
                       
                obj.dxf1 = @(xi)...
                    (0<=xi).*(xi<1).*(...
                        -120*(6*xi.^5 -15*xi.^4 +20*xi.^2 -14)./...
                        (3*xi.^5 -10*xi.^4 +40*xi.^2 +28).^2) +...
                    (1<=xi).*(xi<2).*(...
                          120*(xi-2).^4.*(2*xi+1) ./...
                           (xi.^5 -10*xi.^4 +40*xi.^3 -80*xi.^2 +20*xi -32).^2);    
                
            elseif strcmp(obj.kernel,'Wendland') % very weak for computing Omegaj
                sigma = [3/4; 7/(4*pi); 21/(16*pi)];
                obj.fw = @(q)   sigma(obj.dim)*...
                   (1-q/2).^4 .* (1+2*q);
                obj.fdw =@(q)  sigma(obj.dim)*... 
                  -5*q.*(1-q/2).^3;                  
            else
                error([obj.kernel,' is as kernel not implemented yet']);
            end
            
            
            
            % some preallocation
            obj.Fj_int    = zeros(obj.N,obj.dim);
            obj.Fj_diss   = zeros(obj.N,obj.dim);
            obj.nIj       = zeros(obj.N,obj.dim);
            obj.drhoj_tot = zeros(obj.N,1);
            obj.drhoj_diss = zeros(obj.N,1);
            obj.drhoj_int = zeros(obj.N,1);
            obj.rhoj_half = zeros(obj.N,1);
            obj.vj_half   = zeros(obj.N,obj.dim);
            obj.pj        = zeros(obj.N,1);
            obj.dej       = zeros(obj.N,1);
            obj.ej_half   = zeros(obj.N,1);
            
            obj.pij = [];            
            obj.distances_uptodate = false;
            obj.firststep = true;
            
            obj.cell_of_j      = zeros(obj.N,1);
            obj.obsolete_k_pij = [];


            obj.t = 0;
            %check
            checkDataOnline(obj);
            if obj.N > 50000
                warning('quite a lot particles!')
                keyboard
            end
            obj.IOdata = struct('pjghost',[],... %saved after each mirroring of the particles
                                'vjghost',[],...
                                'ejghost',[],...
                                'rhojghost',[]);
            obj.IOdata_h = struct('pjghost',[],... %saved after each mirroring of the particles
                                'vjghost',[],...
                                'ejghost',[],...
                                'rhojghost',[]);                            
            obj.tmp=0;
            obj.pij_ghost = [];
            obj.Xj_ghost = [];
            obj.compGhost = false; %probably that should be always true (in order to compute volume and pressure)
        end
        %%
        function start_simulation(obj)            
            disp('#### simulation starts ####');
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
                if toc > 3 %seconds
                    disp(['t = ',num2str(obj.t,'%10.5e'),...
                        's (',num2str(round(icount/toc,1)),' iter/sec ,Ncon = ',num2str(size(obj.pij,1)),')']);
                    icount = 1;
                    tic    
                    %disp(['tmp: ',num2str(obj.tmp)]);
                end
                icount = icount +1;                  
                
                %pause
                if any(obj.t > obj.tpause)
                    obj.IO.do(obj,true);  %plot and save data (force it)
                    disp('--- pause (see parameter tpause) ---');
                    keyboard
                    obj.tpause(obj.t>obj.tpause) = inf; %delete this breakpoint
                end
               % obj.checkDataOnline;
            end
            obj.IO.do(obj,true);  %last frame (force it)
            obj.IO.finalize();
            toc(ttic)
        end
        %%        
        function perform_timestep(obj)
             if obj.firststep
                comp_pressure(obj);
                comp_volume(obj); 
             end

             
             % mirror step for boundary condition
             if ~isempty(obj.bc)
                BC(obj)                 
             end

             if obj.compGhost %for data analysis (conservation)
                 obj.IOdata_h.vjghost   = obj.vj_half(obj.Ighost,:);             
                 obj.IOdata_h.ejghost   = obj.ej_half(obj.Ighost);
                 obj.IOdata_h.rhojghost = obj.rhoj_half(obj.Ighost);                               
                 obj.IOdata.vjghost     = obj.vj(obj.Ighost,:);             
                 obj.IOdata.ejghost     = obj.ej(obj.Ighost);
                 obj.IOdata.pjghost     = obj.pj(obj.Ighost,:);             
                 obj.IOdata.rhojghost   = obj.rhoj(obj.Ighost);
             end

             search_neighbours(obj);                                                   
             % mass, momentum and energy fluxes:             
             comp_kernel(obj) 
             if (obj.compOmegaj && ~obj.h_const)
                 comp_Omegaj(obj);                 
             end
             comp_forces(obj)  
             comp_dRho(obj)    
                                        
             if ~obj.isothermal
                comp_de(obj)                             
             end     
             
               
              
             comp_dissipation(obj) %mass, momentum, energy             
               
             %%bc (nrc)
             BCvirtual(obj)  
             
             %delete very small forces
%              epsilon = 1e-14;
%              obj.drhoj(abs(obj.drhoj)<epsilon) = 0;
%              obj.dej(abs(obj.dej)<epsilon) = 0;
%              obj.Fj_tot(abs(obj.Fj_tot)<epsilon) = 0;             
%              
             %experimental:
             %comp_BCdamping(obj)             
             
             update_half_step(obj)
             update_position(obj)   
             update_full_step(obj)
             
             if ~obj.h_const
                 update_h(obj)
             end
             
             comp_pressure(obj);
             comp_volume(obj); %(only for v-scheme, Omegaj, f_ST and phy diss)

        end     
        %%
        function update_dt(obj) 
            if isempty(obj.dt_fix) % only if timestep is not fixed
                if obj.dim == 1
                    obj.dt = obj.dtfactor * min(obj.hj./(obj.cj + abs(obj.vj)));
                else
                    obj.dt = obj.dtfactor * min(obj.hj./(obj.cj + sum(obj.vj(:,1).^2 + obj.vj(:,2).^2,2).^(0.5)));
                end
            end
        end        
        %% neighbour search -  main 
        function search_neighbours(obj)
           if obj.firststep
               disp('create initial cell-structure');               
               initial_search_neighbours(obj);
           elseif (obj.Rtcell < obj.kernel_cutoff * max(obj.hj)) %ToDo: make cells smaller
               obj.Nc = [];  %set cell-division-structure to obsolete
               disp(['create new cell-structure (max(h)= ',num2str( max(obj.hj))]);
               initial_search_neighbours(obj);
           else
               update_neighbours(obj);
           end           
        end        
        %%
        function cell_of_xj = cell_structure(obj,x) % in which cell is x
            %% set cell-division-structure
            if isempty(obj.Nc) %only if old cell-division-structure is obsolete
                 %make cell-division a bit bigger (in order to capture some
                 %derivation in h without creating a new cell-division
                 %structre)
                obj.Rtcell     = obj.kernel_cutoff * max(obj.hj)  * 1.1;
                obj.Nc         = floor((obj.Omega(:,2)-obj.Omega(:,1))/obj.Rtcell);
            end
            % search belonging cell
            NN = size(x,1);
            if obj.dim == 1
                cell_of_xj = floor(ones(NN,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (x - ones(NN,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle      
                if any(cell_of_xj > obj.Nc)  || any(cell_of_xj < 1)                    
                    warning(' - some particles are not in cell structure anymore! - ')
                    keyboard
                end
           
            else
                cell_of_j_2d = floor(ones(NN,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (x - ones(NN,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle
                cell_of_xj = cell_of_j_2d(:,1) + (cell_of_j_2d(:,2)-1).*(obj.Nc(1));    %convert to 1d counting                                                
            end            
            
        end
        %% stop particles on the boundary
        function handleParticlesOnBoundary(obj)
            if obj.dim == 2
                 NN = size(obj.Xj(obj.Icomp),1);
                 cell_of_j_2d = floor(ones(NN,1)*(obj.Nc ./ (obj.Omega(:,2)-obj.Omega(:,1)))'...
                    .* (obj.Xj(obj.Icomp,:) - ones(NN,1)*(obj.Omega(:,1))'))...
                    +1; %in which sector is the particle           
                 %check for particles on the border of Omega:
                on_delOmega = zeros(obj.N,1);
                on_delOmega(obj.Icomp) =...
                       (cell_of_j_2d(:,1) >= obj.Nc(1))...
                     + (cell_of_j_2d(:,1) <= 1 )...
                     + (cell_of_j_2d(:,2) >= obj.Nc(2))...
                     + (cell_of_j_2d(:,2) <= 1) ;
                 
                if any(on_delOmega)                    
%                     warning(' - some particles are not in the cell structure anymore! - ')
%                     keyboard
%                     
                    %or
                    obj.Icomp = obj.Icomp(~ismember(obj.Icomp,find(on_delOmega)));
                    disp('particles on delOmega'); 
%                     keyboard
                end 
            end
        end        
        %%
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
        %%        
        function initial_search_neighbours(obj) %upper right and lower left cells are not included
            %% create cell structure            
            obj.cell_of_j = cell_structure(obj,obj.Xj);
            handleParticlesOnBoundary(obj);
            if obj.dim ==1
                cellshift    = 1;
                Ncell_search = obj.Nc-1;
            else
                cellshift    = [1,obj.Nc(1)-1,obj.Nc(1),obj.Nc(1)+1]; %upper-left side
                Ncell_search = (obj.Nc(1)*obj.Nc(2) - (obj.Nc(1)+1));
            end
            
            if any(obj.cell_of_j > Ncell_search)               
                warning(' - some particles are not in the cell structure! (Increase Omega)- ')
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
        function update_neighbours_new_with_bug (obj) %but a bit faster
 
            %cell-structre:
            cell_of_j_new = cell_structure(obj,obj.Xj);
            handleParticlesOnBoundary(obj);

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
            
           % j_changes_cell = find(cell_shift_of_j); %alternative: a-b ~=0
            

            if obj.dim == 1
                shift = [-1,1,inf,-inf]; 
            else
                shift = [-1,1,...
                       -obj.Nc(1)-1,-obj.Nc(1),-obj.Nc(1)+1,...
                        obj.Nc(1)-1, obj.Nc(1), obj.Nc(1)+1,...
                        inf,-inf]; %todo
            end
            for s = shift %loop over all possible shifts
                j_shift = find(cell_shift_of_j==s);
                new_obsolete_k_pij=[];
                pij_new =[];
                for j = j_shift'    %% loop over all particels with the same shift

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
    
    %other possibility:
    %                     a=logical(obj.AedgesXj(j,:));
    %                     p=obj.pij(a,:);
    %                     m=ismember(p,find(i_in_obs_neighb_cell));
    %                     aa=find(a);
    %                     new_obsolete_k_pij=aa(logical(sum(m,2)))';

                        %find the belonging connectivities                    
                        i_in_obs_neighb_cell(j,1) = 5;  % just a number,to obtain an unsymmetric connection - abs(+-5 -+1)=4
                        temp    = (i_in_obs_neighb_cell'*obj.AedgesXj);                    
                        pij_obs = (abs(temp)==4);
                        NN=size(find(pij_obs'),1);
                        new_obsolete_k_pij(end+(1:NN)) = find(pij_obs');
%                       
%                     else
%                         new_obsolete_k_pij = [];
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
                        pij_new(end+(1:nCon_new),:) = [j*ones(nCon_new,1),i_in_new_neighb_cell]; 
                    else
                        nCon_new = 0;
                    end
                end   
         

                
                if isempty(pij_new)
                    nCon_new = 0;
                else
                    pij_new= pij_new(~ismember(pij_new(:,2),j_shift),:);
                    nCon_new = size(pij_new,1);
                end
                
                if isempty(new_obsolete_k_pij)
                    nCon_ob = 0;
                else
                %check for particles with the same shift (should not be a
                %new neighbour nor an obsolete
              %  keyboard
                   new_obsolete_k_pij = new_obsolete_k_pij (~logical(...
                       sum(ismember(obj.pij(logical(new_obsolete_k_pij),:),j_shift),2)))';                
                   nCon_ob = size(new_obsolete_k_pij,1);
                   
                end
                

                %% apply changes %%
                
                if nCon_ob > 0
                    % clear all obsolete connectivities of node j in the
                    % adjacency matrix nodes <-> edges
                    obj.AedgesXj(:,new_obsolete_k_pij) = 0;
                    %set a self connection for obsoltete particles
                    obj.pij(new_obsolete_k_pij,:) = ...
                                 ones(size(new_obsolete_k_pij,1) ,1) *[1,1]; 
                    obj.obsolete_k_pij = [obj.obsolete_k_pij;...
                                             new_obsolete_k_pij];
                    % http://de.mathworks.com/matlabcentral/newsreader/view_thread/283923
                end
                % insert new connectivies     
                if nCon_new > 0
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
                         obj.pij(end+(1:(nCon_new-nCon_obsolete)),:)= pij_new(nCon_obsolete+1:nCon_new,:);
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
                    
                    for ii = 1:nCon_new  %% improve this!!!
                          obj.AedgesXj (pij_new(ii,2),k_pij_new(ii)) = -1;
                    end
                   % obj.AedgesXj (pij_new(:,2),k_pij_new) = -1*speye(nCon_new); %ingoing  
                    obj.AedgesXj (j , k_pij_new) = 1; %outgoing                          
                end
            end 
                            
            % remove the obsolete ghost particles (which are =inf) in cell_of_j;
            if np_change < 0               
               obj.cell_of_j = obj.cell_of_j(1:end-(-np_change));
               %and shrink connectivity matrix          
               obj.AedgesXj(end-(-np_change)+1:end,:)=[];  
            end
            
%             %check data
%             a=obj.pij;
%             a= a(a(:,1)~=a(:,2),:);
%             initial_search_neighbours(obj);
%             b=obj.pij;
%             if any(any(sort(sort(a,2),1)-sort(sort(b,2),1)))
%                 keyboard
%             end
%            
            
        end                        
        %%  
        function update_neighbours (obj) 
 
            %cell-structre:
            cell_of_j_new = cell_structure(obj,obj.Xj);
            handleParticlesOnBoundary(obj);

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
                 %other possibility:
%                     a=logical(obj.AedgesXj(j,:));
%                     p=obj.pij(a,:);
%                     m=ismember(p,find(i_in_obs_neighb_cell));
%                     aa=find(a);
%                     new_obsolete_k_pij=aa(logical(sum(m,2)))';

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
                         obj.pij(end+(1:(nCon_new-nCon_obsolete)),:)= pij_new(nCon_obsolete+1:nCon_new,:);
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
            
            
            %check data with inital_search algo
%             a=obj.pij;
%             a= a(a(:,1)~=a(:,2),:);
%             initial_search_neighbours(obj);
%             b=obj.pij;
%             if any(any(sort(sort(a,2),1)-sort(sort(b,2),1)))
%                 keyboard
%             end
           
        end                        
        %%         
        function comp_distances(obj)
            obj.xij_h = (obj.Xj(obj.pij(:,1),:)-obj.Xj(obj.pij(:,2),:));
            obj.rij   = sum(obj.xij_h.^2,2).^0.5;
            obj.xij_h = obj.xij_h./(obj.rij*ones(1,obj.dim));  
            
            % flag particle which are within the cutoff radius           ^       
            temp = (obj.rij < obj.kernel_cutoff *...
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
            %% -------------------------------
            function comp_kernel_gauss_hconst(obj)
                sigma = [1/sqrt(pi),1/pi,1/(pi*sqrt(pi))];
                obj.Wij   = 1/(obj.hj(1)^obj.dim)* sigma(obj.dim) *...
                            exp(-(obj.rij(obj.active_k_pij)/obj.hj(1)).^2);
                obj.dWij  = obj.Wij.* -2.*obj.rij(obj.active_k_pij)/obj.hj(1)^2;
                obj.dWij_over_rij = obj.Wij.* -2./obj.hj(1)^2;
                obj.ddWij = -2/obj.hj(1)^2*(obj.dWij.* obj.rij(obj.active_k_pij) + obj.Wij );
            end 
            %% -------------------------------
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
                 
                q = [(obj.rij(obj.active_k_pij)) ./ obj.hj(obj.Ii);
                     (obj.rij(obj.active_k_pij)) ./ obj.hj(obj.Ij)];
 
                Np = sum(obj.active_k_pij);
                I  = q < obj.kernel_cutoff;
                w  = zeros(2*Np,1);            
                w(I) = obj.fw(q(I));

                obj.Wij_hi = w(1:Np)./ (obj.hj(obj.Ii).^obj.dim);
                obj.Wij_hj = w(Np+1:end)./ (obj.hj(obj.Ij).^obj.dim);
                obj.Wij    = 0.5*(obj.Wij_hi + obj.Wij_hj);

                dw = zeros(2*Np,1);
                dw(I) = obj.fdw(q(I,:));

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
                
                
                %axi
                if obj.flag_axi
                    obj.f1j = zeros(obj.N,1);
                    obj.dxf1j = zeros(obj.N,1);

                    I = obj.Xj(:,2) < max(obj.hj)*2;                    
                    obj.f1j = obj.f1(obj.Xj(:,2)./obj.hj(:));
                    obj.dxf1j = obj.dxf1(abs(obj.Xj(:,2)./obj.hj(:)))./obj.hj(:);

                end
             end
        end               
        %%    
        function comp_Omegaj(obj) %very weak with Wendland kernel
            
           hj_bar = 0.5*(obj.hj(obj.Ii)+obj.hj(obj.Ij))  ;
           %normalization:           
           sumW = abs(obj.AedgesXj(:,obj.active_k_pij)) *...
               (obj.Vj(obj.Ii).*obj.Wij)... % fluxes
                +obj.Vj.*obj.fw(0)./ obj.hj.^obj.dim; %himself

           dJ = zeros(size(obj.pij,1),1);
           dJ(obj.active_k_pij,:)=(obj.rhoj(obj.Ii).*obj.Vj(obj.Ii).*obj.Wij); %mj
           Jrhoi = (abs(obj.AedgesXj) * dJ...
                +obj.Vj.*obj.rhoj.*obj.fw(0)./ obj.hj.^obj.dim... %himself...
               )./sumW;
           
           qOmega_dW = zeros(size(obj.pij,1),1);
           qOmega_W  = zeros(size(obj.pij,1),1);

           q = (obj.rij(obj.active_k_pij)) ./ hj_bar;
           dWdhij_dW  = - q.*obj.dWij;
           dWdhij_W   = - obj.dim./hj_bar .*obj.Wij; 
           
           dhdrho =  -1/obj.dim * obj.eta * obj.mj.^(1/obj.dim).*...
                            obj.rhoj.^(-1/obj.dim -1);  
           
           if strcmp(obj.scheme,'m')                  
               qOmega_dW(obj.active_k_pij,:)  = obj.Vj(obj.Ij).*dWdhij_dW.*(obj.rhoj(obj.Ij)-Jrhoi(obj.Ij));
               qOmega_W(obj.active_k_pij,:)   = obj.Vj(obj.Ij).*dWdhij_W.*(obj.rhoj(obj.Ij)-Jrhoi(obj.Ij));
               

               J1 = ((abs(obj.AedgesXj) * qOmega_W)... %add up all fluxes (no distinction between in and outcomming)->abs                                      
                            - obj.Vj.*(obj.rhoj-Jrhoi).*obj.dim*obj.fw(0)./ obj.hj.^(obj.dim+1)...  % consider selfinteraction;
                            );                                 
               J2 = (abs(obj.AedgesXj) * qOmega_dW);
              
               obj.Omegaj = 1  - (dhdrho.*(J1+J2)./sumW);                            


           elseif strcmp(obj.scheme,'v')
               warning('no normalization!')
               qOmega(obj.active_k_pij,:) = obj.Vj(obj.Ij).*dWdhij_hi;

               obj.Omegaj = 1  - (obj.rhoj.*dhdrho.*...
                            (...
                            (abs(obj.AedgesXj) * qOmega)... %add up all fluxes (no distinction between in and outcomming)->abs
                            - obj.Vj.*obj.dim*obj.fw(0)./ obj.hj.^(obj.dim+1)...  % consider selfinteraction;
                            ))./sumW;   
           elseif strcmp(obj.scheme,'n')
               warning('no normalization!')
               qOmega(obj.active_k_pij,:) = dWdhij_hi;

               obj.Omegaj = 1  - (obj.mj.*dhdrho.*...
                            (...
                            (abs(obj.AedgesXj) * qOmega)... %add up all fluxes (no distinction between in and outcomming)->abs
                            - obj.dim*obj.fw(0)./ obj.hj.^(obj.dim+1)...  % consider selfinteraction;
                            ))./sumW;                 
           else
               error('only m,n,v-schemes are supported right now');
           end
 
%            disp(['min:',num2str(min(obj.Omegaj)),' max:',num2str(max(obj.Omegaj))]);
%            if min(obj.Omegaj)<0
%                keyboard
%            end
        end
        %%
        function update_full_step(obj)             
            if ~obj.firststep
                obj.rhoj(obj.Icomp) = obj.rhoj_half(obj.Icomp)+ 0.5*obj.dt * obj.drhoj_tot(obj.Icomp);
                obj.vj(obj.Icomp,:)  = obj.vj_half(obj.Icomp,:) + 0.5*obj.dt *...
                    obj.Fj_tot(obj.Icomp,:)./(obj.mj(obj.Icomp)*ones(1,obj.dim));
                if ~isempty(obj.g_ext)
                    obj.vj(obj.Icomp,:)  = obj.vj_half(obj.Icomp,:) + 0.5*obj.dt * ones(size(obj.Icomp))*obj.g_ext;
                end
                if ~obj.isothermal
                      obj.ej(obj.Icomp) = obj.ej_half(obj.Icomp)+ 0.5*obj.dt * obj.dej(obj.Icomp);
                end                
            end
        end
        %%
        function update_half_step(obj) 
            if obj.firststep
                obj.rhoj_half(obj.Icomp) = obj.rhoj(obj.Icomp)  + 0.5*obj.dt*obj.drhoj_tot(obj.Icomp);
                obj.vj_half(obj.Icomp,:)  = obj.vj(obj.Icomp,:) + 0.5*obj.dt * ...
                    (obj.Fj_tot(obj.Icomp,:))./(obj.mj(obj.Icomp)*ones(1,obj.dim));
                if ~isempty(obj.g_ext)
                    obj.vj_half(obj.Icomp,:)  = obj.vj(obj.Icomp,:) + 0.5*obj.dt * ones(size(obj.Icomp))*obj.g_ext;
                end
                if ~obj.isothermal
                    obj.ej_half(obj.Icomp) = obj.ej(obj.Icomp)  + 0.5*obj.dt*obj.dej(obj.Icomp);
                end
                
                obj.firststep=false;
            else               
                obj.rhoj_half(obj.Icomp) = obj.rhoj_half(obj.Icomp)+ obj.dt * obj.drhoj_tot(obj.Icomp);
                obj.vj_half(obj.Icomp,:)  = obj.vj_half(obj.Icomp,:) + obj.dt *...
                    (obj.Fj_tot(obj.Icomp,:))./(obj.mj(obj.Icomp)*ones(1,obj.dim));
                if ~isempty(obj.g_ext)
                    obj.vj_half(obj.Icomp,:)  = obj.vj_half(obj.Icomp,:) + obj.dt * ones(size(obj.Icomp))*obj.g_ext;
                end
                if ~obj.isothermal
                     obj.ej_half(obj.Icomp) = obj.ej_half(obj.Icomp)+ obj.dt * obj.dej(obj.Icomp);
                end
            end
            
        end
        %%
        function update_position(obj) 
             obj.Xj(obj.Icomp,:)      = obj.Xj(obj.Icomp,:) + obj.dt * obj.vj_half(obj.Icomp,:);
             obj.distances_uptodate = false;
        end
        %%
        function update_h(obj)
           obj.hj(obj.Icomp) = obj.hj(obj.Icomp)...
               - 1/obj.dim * obj.dt * (obj.hj(obj.Icomp)./obj.rhoj(obj.Icomp)) .* obj.drhoj_tot(obj.Icomp);
           %disp(['min: ',num2str(min(obj.hj)),'; max: ',num2str(max(obj.hj))]);
        end
        %%
        function comp_volume(obj)
            %only necessary for v-scheme
            obj.Vj = obj.mj(:)./obj.rhoj(:);
        end        
        %%
        function comp_pressure(obj)
            %% -----------------------
            function comp_pressure_EOS_isothermal(obj)
                obj.pj(obj.Icomp,:) = obj.c0j(obj.Icomp) .* obj.c0j(obj.Icomp) ...
                                     .*(obj.rhoj_real(obj.Icomp) - obj.rho0j(obj.Icomp));
                %c= sqrt(dp/drho)=c0 -> no change
            end
            %%
            function comp_pressure_EOS_mie_gruneinsen(obj) %double check that
                a0 = obj.rho0j(obj.Icomp).*obj.c0j(obj.Icomp).^2; 
                b0 = a0.*(1+2*(obj.MG_Sj(obj.Icomp)-1));
                c0 = a0.*(2*(obj.MG_Sj(obj.Icomp)-1)+3*(obj.MG_Sj(obj.Icomp)-1).^2);

                etaMG = obj.rhoj_real(obj.Icomp)./obj.rho0j(obj.Icomp)-1;

                pH = a0.*etaMG + (etaMG>0) .* (b0.*etaMG.^2 + c0.*etaMG.^3);

                obj.pj(obj.Icomp,:) = (1 - 0.5*obj.MG_Gammaj(obj.Icomp).*etaMG) .*pH ...
                                    + obj.MG_Gammaj(obj.Icomp) .*obj.rhoj_real(obj.Icomp).* obj.ej(obj.Icomp); 
                
                % update speed of sound: (copy paste from LimeSPH)
                obj.cj(obj.Icomp,:) = (-0.5*pH./obj.rho0j(obj.Icomp,:) ...
                                    + (1-0.5*obj.MG_Gammaj(obj.Icomp).*etaMG) .*...
                                    (a0 + ...
                                    (2*b0.*etaMG + 3*c0.*etaMG.^2).*(etaMG>0)...                                          
                                     )./obj.rho0j(obj.Icomp,:) +...
                                     obj.MG_Gammaj(obj.Icomp).*obj.ej(obj.Icomp,:)).^0.5; 
            end
            %%
            function comp_pressure_EOS_water(obj)
                obj.pj(obj.Icomp,:) = obj.c0j(obj.Icomp) .* obj.c0j(obj.Icomp).* ...
                                    obj.rho0j(obj.Icomp).*...
                                     ((obj.rhoj_real(obj.Icomp)./obj.rho0j(obj.Icomp)).^7 - 1)/7;
            end
            %%
            function comp_pressure_EOS_IdealGas(obj,Gamma)
                obj.pj(obj.Icomp,:) = (Gamma-1)*obj.rhoj_real(obj.Icomp,:).*obj.ej(obj.Icomp,:);                
                % from paper %c=sqrt(Gamma*p/q)
                 obj.cj(obj.Icomp,:) = sqrt(Gamma*obj.pj(obj.Icomp)./obj.rhoj_real(obj.Icomp));
               
%                obj.cj(obj.Icomp,:) = sqrt((gamma-1)*obj.ej(obj.Icomp));
            end
            
            
            %% -----------------------
            %compute 3d-velocity in the case of an axisymetric simulation
            if obj.flag_axi
                obj.rhoj_real = obj.rhoj./(2*pi*obj.Xj(:,2));
            else
                obj.rhoj_real = obj.rhoj;
            end
            
            if strcmp(obj.EOS,'ISO')
                comp_pressure_EOS_isothermal(obj)
            elseif strcmp(obj.EOS,'MG')
                comp_pressure_EOS_mie_gruneinsen(obj)
            elseif strcmp(obj.EOS,'Water')
                comp_pressure_EOS_water(obj)
            elseif strcmp(obj.EOS,'IdealGas53')
                Gamma53= 5/3;
                comp_pressure_EOS_IdealGas(obj,Gamma53)
            elseif strcmp(obj.EOS,'IdealGas14')
                Gamma53= 1.4;
                comp_pressure_EOS_IdealGas(obj,Gamma53)
            else
                error('EOS not supported');
            end
        end
        %% change of energy (ideal process -reversible and adiabatic)
        function comp_de(obj)
            obj.dej = 0*obj.dej ;
            if ~obj.flag_axi
                obj.dej(obj.Icomp) = obj.pj(obj.Icomp)./...
                             ((obj.rhoj(obj.Icomp)).^2) .* obj.drhoj_int(obj.Icomp);
            else
               obj.dej(obj.Icomp) =2*pi*obj.pj(obj.Icomp)./ obj.rhoj(obj.Icomp) .*...
                             (-obj.vj(obj.Icomp,2) + ...
                              obj.Xj(obj.Icomp,2)./obj.rhoj(obj.Icomp) .* obj.drhoj_int(obj.Icomp));                                   
            end
        end          
        %%
        function comp_dRho(obj)
            %%--------------------
            function comp_dRho_m_scheme(obj)
                obj.drhoj_int  = 0*obj.drhoj_int; 
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
                obj.drhoj_int = ((obj.AedgesXj>0) * qrho_ij +...
                            ( obj.AedgesXj<0) * qrho_ji)./obj.Omegaj;  
            end
            %%
            function comp_dRho_n_scheme(obj)
                obj.drhoj_int  = 0*obj.drhoj_int; 
                qrho_ij = zeros(size(obj.pij,1),1);
                qrho_ji = zeros(size(obj.pij,1),1);
                
                qrho_ij(obj.active_k_pij,:) =  ... % i-> j %hi
                    sum(... %scalar product v*n 
                    (obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)) .*...
                    (obj.dWij_hi*ones(1,obj.dim)).*obj.xij_h((obj.active_k_pij),:)...
                    ,2);
                qrho_ji(obj.active_k_pij,:) =  ... % j -> i %hj
                     sum(... %scalar product v*n 
                    (obj.vj(obj.Ij,:)-obj.vj(obj.Ii,:)) .*...
                    (obj.dWij_hj*ones(1,obj.dim)).*-obj.xij_h((obj.active_k_pij),:)...
                    ,2);
                
                
                %add up all the corresponding density flux in each node
                obj.drhoj_int = ((obj.AedgesXj>0) * qrho_ij +...
                            ( obj.AedgesXj<0) * qrho_ji).*obj.mj./obj.Omegaj;  
            end
            %%
            function comp_dRho_v_scheme(obj)
                obj.drhoj_int  = 0*obj.drhoj_int; 
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
                obj.drhoj_int = ((obj.AedgesXj>0) * qrho_ij +...
                             (obj.AedgesXj<0) * qrho_ji)...
                            .*obj.rhoj./obj.Omegaj;
            end
            %%
            function comp_dRho_axi(obj)
                obj.drhoj_int  = 0*obj.drhoj_int; 
                qrho_ij = zeros(size(obj.pij,1),1);
                qrho_ji = zeros(size(obj.pij,1),1);
                                
                
                %first term (r-terms)
                qrho_ij(obj.active_k_pij,1) =  ... % i-> j %hi
                    obj.mj(obj.Ij) .*...                     
                    (obj.f1j(obj.Ii).*obj.vj(obj.Ii,2) - obj.f1j(obj.Ij).*obj.vj(obj.Ij,2)) .*...
                    (obj.dWij_hi).*obj.xij_h((obj.active_k_pij),2);                    
                
                
                qrho_ji(obj.active_k_pij,:) =  ... % j -> i %hj
                    obj.mj(obj.Ii).*...
                    (obj.f1j(obj.Ii).*obj.vj(obj.Ij,2) - obj.f1j(obj.Ij).*obj.vj(obj.Ii,2)) .*...
                    obj.dWij_hj.*-obj.xij_h((obj.active_k_pij),2);
                
                
                %add up all the corresponding density fluxes in each node
                obj.drhoj_int = ((obj.AedgesXj>0) * qrho_ij +...
                            ( obj.AedgesXj<0) * qrho_ji);                 
                
                % second term (still r-terms, but now symmetric)
                qrho_ij(obj.active_k_pij,1) = ...
                    obj.mj(obj.Ij).*...
                    (obj.dxf1j(obj.Ii).*obj.vj(obj.Ii,2) - obj.dxf1j(obj.Ij).*obj.vj(obj.Ij,2)).*...
                    obj.Wij;
                obj.drhoj_int = obj.drhoj_int +...
                    obj.AedgesXj * qrho_ij;
                        
                %third term (z-term)
                qrho_ij(obj.active_k_pij,1) =  ... % i-> j %hi
                    obj.mj(obj.Ij) .*...                     
                    (obj.vj(obj.Ii,1) - obj.vj(obj.Ij,1)) .*...
                    (obj.dWij_hi).*obj.xij_h((obj.active_k_pij),1);                    
                
                
                qrho_ji(obj.active_k_pij,:) =  ... % j -> i %hj
                    obj.mj(obj.Ii).*...
                    (obj.vj(obj.Ij,2) - obj.vj(obj.Ii,2)) .*...
                    obj.dWij_hj.*-obj.xij_h((obj.active_k_pij),2);
              
                obj.drhoj_int = obj.drhoj_int + obj.f1j.*...
                    ((obj.AedgesXj>0) * qrho_ij +...
                            ( obj.AedgesXj<0) * qrho_ji);                             
            end
            
            
            %% ------------------
                       
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            if strcmp(obj.scheme,'m')            
                comp_dRho_m_scheme(obj);
            elseif strcmp(obj.scheme,'n')
                comp_dRho_n_scheme(obj);
            elseif strcmp(obj.scheme,'v')
                comp_dRho_v_scheme(obj);
            elseif strcmp(obj.scheme,'a')
                comp_dRho_axi(obj);
            else
                error('scheme not supported');
            end    
            obj.drhoj_tot = obj.drhoj_int;
        end     
        %%
        function comp_forces(obj)
                    
            %%--------------------
            function comp_Fint_m_scheme(obj) %internal forces  
                %compute the flux
                qFij_int = zeros(size(obj.pij,1),1);
                p_rho_omega_j = obj.pj./(obj.Omegaj.*obj.rhoj.*obj.rhoj);
                
                qFij_int(obj.active_k_pij,:)=...
                    -(obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... % hi
                    (p_rho_omega_j(obj.Ii).*obj.dWij_hi)...  %<- change here for variable h
                    +...
                    obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... %hj
                    (p_rho_omega_j(obj.Ij).*obj.dWij_hj)...  %<- change here
                    );      
                
                 %spread fluxes to the nodes       
                 obj.Fj_int = obj.AedgesXj * ((qFij_int*ones(1,obj.dim)).*obj.xij_h);
            end            
            %%
            function comp_Fint_n_scheme(obj) %internal forces  
                %compute the flux
                qFij_int = zeros(size(obj.pij,1),1);
                pVV_omega_j = obj.Vj.^2.*obj.pj./(obj.Omegaj);

                qFij_int(obj.active_k_pij,:)=...
                    -(...
                    (pVV_omega_j(obj.Ii).*obj.dWij_hi)... 
                    +...
                    (pVV_omega_j(obj.Ij).*obj.dWij_hj)...  
                    );           
                  %spread fluxes to the nodes       
                 obj.Fj_int = obj.AedgesXj * ((qFij_int*ones(1,obj.dim)).*obj.xij_h);
            end
            %%
            function comp_Fint_v_scheme(obj) %internal forces  
                %compute the flux
                qFij_int = zeros(size(obj.pij,1),1);
                p_omega_j = obj.pj./(obj.Omegaj);

                qFij_int(obj.active_k_pij,:)=...
                    -(obj.Vj(obj.Ii).*obj.Vj(obj.Ij) .*... % hi
                    (p_omega_j(obj.Ii).*obj.dWij_hi)...  
                    +...
                    obj.Vj(obj.Ii).*obj.Vj(obj.Ij) .*... %hj
                    (p_omega_j(obj.Ij).*obj.dWij_hj)... 
                    );
                  %spread fluxes to the nodes       
                 obj.Fj_int = obj.AedgesXj * ((qFij_int*ones(1,obj.dim)).*obj.xij_h);

            end
            %%
            function comp_Fj_phy_diss(obj)      %page415-violeau | particle-friction   
                %compute the flux
                warning('check this')
                qFj_ph_diss_ij = zeros(size(obj.pij,1),1);
                qFj_ph_diss_ij(obj.active_k_pij,:)  = 2*(obj.dim+2)*obj.mu*(...
                    obj.Vj(obj.Ii,:) .* obj.Vj(obj.Ij,:).*... Vi*Vj
                             sum(...
                             ((obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)).*... $(vi-vj)*hat(xij)
                             obj.xij_h(obj.active_k_pij,:))...
                             ,2).*...
                             obj.dWij);  %dWij
                 %spread fluxes to the nodes       
                 obj.Fj_phy_diss = obj.AedgesXj * ((qFj_ph_diss_ij*ones(1,obj.dim)).*obj.xij_h);                            
            end
            %%
            function comp_Fj_ST(obj)  %Violeau423        
                 comp_normal(obj)
                 warning('check this')

                %compute the flux
                 qFj_ST_ij = zeros(size(obj.pij,1),obj.dim);
                 qFj_ST_ij(obj.active_k_pij,:)  =  - obj.beta* ((obj.Vj(obj.Ii,:) .* obj.Vj(obj.Ij,:)) ... Vi*Vj
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
                obj.Fj_ST  = obj.AedgesXj * qFj_ST_ij;         
            end
            
            function comp_Fj_axi(obj) %dWij_hi /hj ?
                %compute the flux
                qFij_int = zeros(size(obj.pij,1),1);
                %symmetry is at y=0
                p_r_rho = obj.pj.*obj.Xj(:,2)./(obj.rhoj.*obj.rhoj).*obj.f1j;
                
                qFij_int(obj.active_k_pij,:)= -2*pi*...
                    (obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... % hi
                    (p_r_rho(obj.Ii).*obj.dWij_hi)...  %<- change here for variable h
                    +...
                    obj.mj(obj.Ii).*obj.mj(obj.Ij) .*... %hj
                    (p_r_rho(obj.Ij).*obj.dWij_hj)...  %<- change here
                    );      
                                               
                
                 %spread fluxes to the nodes       
                 obj.Fj_int = obj.AedgesXj * ((qFij_int*ones(1,obj.dim)).*obj.xij_h);                
                 
                 %add hoop stress and correction term
                  obj.Fj_int(:,2) = obj.Fj_int(:,2) + ...
                      (2*pi*obj.pj(:)./obj.rhoj(:)).* ...
                      (1- obj.Xj(:,2)./obj.f1j .* obj.dxf1j);
            end
            
            %%--------------------
            
            if ~obj.distances_uptodate
                comp_distances(obj);
            end
            %internal forces
            if strcmp(obj.scheme,'m')            
                comp_Fint_m_scheme(obj);
            elseif strcmp(obj.scheme,'n')
                comp_Fint_n_scheme(obj);
            elseif strcmp(obj.scheme,'v')
                comp_Fint_v_scheme(obj);
            elseif strcmp(obj.scheme,'a')
                comp_Fj_axi(obj);
            else
                error('selected scheme is not supported');
            end
            obj.Fj_tot = obj.Fj_int;            
            %physical dissipation
            if obj.mu ~=0
                comp_Fj_phy_diss(obj);
                obj.Fj_tot(obj.Icomp,:) = obj.Fj_tot(obj.Icomp,:) + obj.Fj_phy_diss(obj.Icomp,:);
            end
            
            % surface tension
            if obj.beta ~= 0
                comp_Fj_ST(obj)
                obj.Fj_tot(obj.Icomp,:) = obj.Fj_tot(obj.Icomp,:) + obj.Fj_ST(obj.Icomp,:);
            end
        end
        %%
        function comp_dissipation(obj) %artificial
            % parameters:
            alpha_mass      = obj.art_diss_para.alpha_mass;
            beta_mass       = obj.art_diss_para.beta_mass;
            alpha_viscosity = obj.art_diss_para.alpha_viscosity;
            beta_viscosity  = obj.art_diss_para.beta_viscosity;
            alpha_energy    = obj.art_diss_para.alpha_energy;           
            beta_energy     = obj.art_diss_para.beta_energy;
            % reset data:
            obj.drhoj_diss = 0*obj.drhoj_diss; 
            qrho_ij_diss   = zeros(size(obj.pij,1),1);
            qFij_diss      = zeros(size(obj.pij,1),1);           
            obj.dej_diss   = 0*obj.dej_diss; 
            qe_ij_diss     = zeros(size(obj.pij,1),1);
            
            %% some general quantities:
            vijxijh = sum(...
                     (obj.vj(obj.Ii,:)-obj.vj(obj.Ij,:)).*... (vi-vj)*hat(xij)
                     obj.xij_h(obj.active_k_pij,:)...
                     ,2);
            rho_ij_bar = (obj.rhoj(obj.Ii) + obj.rhoj(obj.Ij))/2;
            cij_bar    = (obj.cj(obj.Ii)   + obj.cj(obj.Ij)  )/2;

            %% Dissipative mass flux (Zizis2014)
            v_sig_rho = (cij_bar...
                   - 0.5*beta_mass*vijxijh);
            v_sig_rho(vijxijh>0)=inf; %consider only aproaching particles

                            
            % compute fluxes
            qrho_ij_diss(obj.active_k_pij,:) = alpha_mass* ... % j->i
                obj.mj(obj.Ii).* obj.mj(obj.Ij).*...
                (obj.pj(obj.Ii)-obj.pj(obj.Ij))./...
                (rho_ij_bar .* v_sig_rho)...
                .*obj.dWij;
            %note: hat(xij)*gradiWij = hat(xij)*Wij'*hat(xij)) = Wij' =
            %hat(xij)*gradjWji
          
            
            % spread flux onto nodes
            obj.drhoj_diss = (obj.AedgesXj*qrho_ij_diss)./obj.mj;

            obj.drhoj_tot = obj.drhoj_tot + obj.drhoj_diss;

            
            %% viscosity - dissipative velocity - Iason              
             v_sig_force = (cij_bar...
                    - 0.5*beta_viscosity*vijxijh);                   
             v_sig_force(vijxijh>0) = 0;
             %compute the fluxes
             qFij_diss(obj.active_k_pij,:)  = ...
                 ((obj.mj(obj.Ii) .* obj.mj(obj.Ij).* ...
                 alpha_viscosity .*v_sig_force .* vijxijh) ./...
                 rho_ij_bar...
                 .*obj.dWij);
             
             % spread flux onto nodes
             obj.Fj_diss = obj.AedgesXj * ((qFij_diss*ones(1,obj.dim)).*obj.xij_h);                        
             obj.Fj_tot = obj.Fj_tot + obj.Fj_diss;                        
            if ~obj.isothermal % if energy dependent
                %--------------------            
                %energy dissipation because of the viscosity:
                de = - sum(obj.Fj_diss .*obj.vj,2)./obj.mj;
                obj.dej(obj.Icomp,:) = obj.dej(obj.Icomp,:)+ de(obj.Icomp,:);
                %--------------
                
                %% thermal conductivity (Zizis2014)

                %approach 1           
                % zizis:
                v_sig_energy = cij_bar - 0.5*beta_energy*vijxijh;   
                % price:
%                 v_sig_energy = (abs(obj.pj(obj.Ii)-obj.pj(obj.Ij))./rho_ij_bar).^0.5;
                
                v_sig_energy(vijxijh>0)=0;
                
                % compute dissipative fluxes :
                qe_ij_diss(obj.active_k_pij,:) = alpha_energy* ... % j-> i
                    v_sig_energy.*...
                    obj.mj(obj.Ii).*obj.mj(obj.Ij).*...
                    (obj.ej(obj.Ii)-obj.ej(obj.Ij))./...
                    (rho_ij_bar )...
                    .*obj.dWij;
                
                % spread fluxes onto nodes
                obj.dej_diss =  (obj.AedgesXj *  qe_ij_diss)./obj.mj;
  
                obj.dej(obj.Icomp,:)   = obj.dej(obj.Icomp,:) +...
                                             obj.dej_diss(obj.Icomp,:);
            end
        end
        
        %% remove?
        function comp_dc4j(obj)
            % compute half-step vaules first
            comp_characteristics_half(obj);
            
            
            obj.dc4j  = 0*obj.dc4j; 
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
            obj.drhoj_int = ((obj.AedgesXj>0) * qrho_ij +...
                        ( obj.AedgesXj<0) * qrho_ji)./obj.Omegaj;  

        end
        %% remove?
        function comp_characteristics_half(obj)           
            % one should use pj_half!!!
            I=obj.Iin;
            
            rho_ref = 1;
            p_ref = 0;
            c_ref = 1;
            v_ref = 0;
            n = 1;
            
            %vertical perturbation
            dv = (obj.vj(I,:) - v_ref) - (obj.vj(I,:) - v_ref)*n';
            
            %with transformation:
            obj.c1j_half= -c_ref.^2 .*(obj.rhoj_half(I) - rho_ref)...
                + (obj.pj(I) - p_ref);
            obj.c2j_half= obj.rhoj(I).*c_ref .* dv;
            obj.c3j_half= obj.rhoj(I).*c_ref .* (obj.vj(I,:) - v_ref)*n'...
                + (obj.pj(I) - p_ref);
%             obj.c4j= - obj.rhoj(I).*c_ref .* ((obj.vj(I,:) - v_ref)*n')...
%             + (obj.pj(I) - p_ref);
            
            % per integration
            obj.c4j_half = obj.c4j(I) + 0.5 * obj.dc4j(I);
        end        
        %%
        function BC(obj)
            resetData(obj)
            for iboun = 1:size(obj.bc,2)
                boun = obj.bc(iboun);
                if strcmp (boun.type,'nrc') %characteristic approch
                    BCnrc_init (obj,iboun) 
                elseif strcmp (boun.type,'nrp') %characteristic approch with particles
                    BCnrc_init (obj,iboun)  %same init-function
                    BCnrp (obj, iboun);
                elseif strcmp (boun.type,'nrm') %mirror approch
                    BCnrm (obj,iboun);
                elseif strcmp (boun.type,'noflow') 
                    BCnoflow(obj,iboun);
                elseif strcmp (boun.type,'cut')
                    %do nothing
                else                     
                    error('no such boundary condition implemented');
                end
            end                  
        end
        %%          
        function BCnrc_init (obj,iboun)
            if obj.firststep %search for the boundary particles and set reference (one time in the beginning)
               dx = mean(obj.Vj)^(1/obj.dim);
%               obj.bc(iboun).kb = obj.point_to_line(obj.Xj(obj.Iin,:),obj.bc(iboun).p1,obj.bc(iboun).p2)< 3*obj.eta*dx;
               n =  obj.bc(iboun).outer_normal;
               bp = obj.bc(iboun).bp;
               obj.bc(iboun).kb = abs(sum((obj.Xj- (ones(obj.N,1)*bp)) .* ...
                                             (ones(obj.N,1)*n) ...
                                         ,2)) < 3*obj.eta*dx;
               obj.bc(iboun).rho_ref = obj.rhoj(obj.bc(iboun).kb);
               obj.bc(iboun).p_ref   = obj.pj(obj.bc(iboun).kb);
               obj.bc(iboun).v_ref   = obj.vj(obj.bc(iboun).kb,:); 
%                obj.bc(iboun).v_ref  = obj.bc(iboun).v_ref *0;
            end
        end    
        %%
        function BCvirtual(obj)
            for i = 1:size(obj.bc,2)
                boun = obj.bc(i);
                if strcmp (boun.type,'nrc') %characteristic approach
                    BCnrc (obj,i);
                end
            end                  
        end           
        %%
        function BCnrc (obj,iboun)
            kbb = (find (obj.bc(iboun).kb));

            qdW = zeros(size(obj.pij,1),obj.dim);
            qdW(obj.active_k_pij,:)= (obj.mj(obj.Ii).*...
                        obj.dWij_hi*...
                        ones(1,obj.dim)).*...
                        obj.xij_h(obj.active_k_pij,:);

            sumdW = obj.AedgesXj(kbb,:) * qdW;
            
            v_ref = obj.bc(iboun).v_ref;
            p_ref = obj.bc(iboun).p_ref;
            rho_ref = obj.bc(iboun).rho_ref;
            c_ref = obj.c0j(kbb,:);

            %or dW normal?
            n = obj.bc(iboun).outer_normal;
          %  nW = -sumdW(kbb,:)./((sum(sumdW(kbb,:).^2,2).^0.5)*ones(1,obj.dim));
            n_len = sum(sumdW .*(ones(size(kbb,1),1)*n),2);           
            
            % ----------------------------------------------
%             %smoothend values:                                                
%             NN=size(kbb,1);
%            for k = 1:NN   %evaluate point for point             %Violeau (5.86)
%                 II1 = (ismember(obj.pij(obj.active_k_pij,1),kbb(k)));
%                 II2 = (ismember(obj.pij(obj.active_k_pij,2),kbb(k)));
%                 if (sum(II1)+sum(II2))==0
%                     continue
%                 end
%                 %normalization constant
%                 sumW = sum(obj.Vj(obj.pij(II1,2)).*obj.Wij(II1),1)...
%                       +sum(obj.Vj(obj.pij(II2,1)).*obj.Wij(II2),1);
%                 if sumW == 0 % no data to sum up, so just let it be zero
%                     continue              
%                 end
%                 dat =obj.rhoj;
%                 dat=dat.*obj.Vj;
%                 rhojkbb(k,:) = (sum((obj.Wij(II1)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II1,2),:)),...
%                                     1)+...
%                                 sum((obj.Wij(II2)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II2,1),:)),...
%                                     1))...
%                                     /sumW; %sum all rows
%                 dat =obj.pj;
%                 dat=dat.*obj.Vj;
%                 pjkbb(k,:) = (sum((obj.Wij(II1)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II1,2),:)),...
%                                     1)+...
%                                 sum((obj.Wij(II2)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II2,1),:)),...
%                                     1))...
%                                     /sumW; %sum all rows
%                 dat =obj.vj;
%                 dat=dat.*(obj.Vj*ones(1,obj.dim));
%                 vjkbb(k,:) = (sum((obj.Wij(II1)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II1,2),:)),...
%                                     1)+...
%                                 sum((obj.Wij(II2)* ones(1,size(dat,2))...
%                                     .* dat(obj.pij(II2,1),:)),...
%                                     1))...
%                                     /sumW; %sum all rows
%            end
           % ----------------------------------------------
%             keyboard
            %point values
            rhojkbb = obj.rhoj(kbb);
            vjkbb = obj.vj(kbb,:);
            pjkbb=obj.pj(kbb);
            
            
            J1 = -c_ref.^2 .*(rhojkbb - rho_ref)...
                + (pjkbb - p_ref);
            %J1 very small
            J2 = obj.rhoj(kbb).*c_ref .* ((vjkbb - v_ref)*n')...
                + (pjkbb - p_ref);

            J3 = - obj.rhoj(kbb).*c_ref .* ((vjkbb - v_ref)*n')...
            + (pjkbb - p_ref);
%             if any(J2)
%                 y=J1';
%                 x=obj.Xj(kbb)';
%                 xx=obj.Xj(end)+(obj.Xj(end)-obj.Xj(end-1));
%                 interp1(x,y,xx)
%                 keyboard
%             end

%             J2 = 2*J2(end)-J2(end-1);
%             J1 = 2*J1(end)-J1(end-1);

%             J2=mean(J2);
%             J1=mean(J1);

             % possibility 1: n= n
            J3=J3*0;         
            % possibility 2: -> n= -n
%                     J1=J1*0;
%                     J2=J2*0;    


           % virtual values
            rhojb = rho_ref + ...
                1./(c_ref.^2).*(-J1+0.5*J2+0.5*J3);
            vjb = v_ref + ...
                (1./(2*c_ref.*rhojb) .* (J2 - J3))*n;
            pjb = p_ref+...
                    0.5*(J2+J3);

                %tangential:            
%             vjb = vjb + obj.vj(kbb,:) - (obj.vj(kbb,:)*n')*n;
                
            if obj.compOmegaj
                omega = obj.Omegaj(kbb);
            else
                omega = 1;
            end
            
            % dissipative mass flux
            obj.bc(iboun).drhoj_diss = sum((obj.vj(kbb,:)-vjb) .* (n_len*n),2)./omega;
            obj.bc(iboun).drhoj_diss = obj.bc(iboun).drhoj_diss;
            
            obj.drhoj_tot(kbb) = obj.drhoj_tot(kbb) - obj.bc(iboun).drhoj_diss;
            
            % dissipative force
            obj.bc(iboun).dvj_diss = -(((obj.pj(kbb)./obj.rhoj(kbb).^2 + pjb./rhojb.^2./omega)...
                *ones(1,obj.dim))...
                .*(n_len*n));
            obj.bc(iboun).dvj_diss = obj.bc(iboun).dvj_diss;
            
            
            obj.Fj_tot(kbb,:) = obj.Fj_tot(kbb,:) - obj.bc(iboun).dvj_diss .* (obj.mj(kbb)*ones(1,obj.dim));
            

            
            % dissipative energy
            obj.bc(iboun).dej_diss = obj.pj(kbb)./obj.rhoj(kbb).^2 .* obj.bc(iboun).drhoj_diss;
            obj.dej(kbb) = obj.dej(kbb) -  obj.bc(iboun).dej_diss;
        end
        %%
        function BCnrm (obj,iboun)
            %%
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
            %%
            
            if obj.firststep %search for the boundary particles

                dx = mean(obj.Vj)^(1/obj.dim);                
%                 obj.bc(iboun).kb = obj.point_to_line(obj.Xj(obj.Iin,:),obj.bc(iboun).p1,obj.bc(iboun).p2)< dx;
                n  = obj.bc(iboun).outer_normal;
                bp = obj.bc(iboun).bp;
                NN = size(obj.Iin,1);
                obj.bc(iboun).kb = abs(sum((obj.Xj(obj.Iin,:)- (ones(NN,1)*bp)) .* (ones(NN,1)*n),2)) < dx;

                % temp
               obj.bc(iboun).rho_ref = obj.rhoj(obj.bc(iboun).kb);
               obj.bc(iboun).p_ref   = obj.pj(obj.bc(iboun).kb);
               obj.bc(iboun).v_ref   = obj.vj(obj.bc(iboun).kb,:); 
               obj.bc(iboun).e_ref   = obj.ej(obj.bc(iboun).kb); 

            end
            
            kb = obj.bc(iboun).kb; %boundary points
            ki = find(kb==0);      %inner points

            outer_normal = obj.bc(iboun).outer_normal;

            %change of sign: (set this values to 1 for second
            %approach - only mirroring the particles)
            if obj.dim == 1
                v_factor = -1;
            else
                v_factor = [-1,-1];
            end
            p_factor = -1;
            v0   =  obj.bc(iboun).v_ref(1,:);
            e0   =  obj.bc(iboun).e_ref(1,:);
            p0   =  obj.bc(iboun).p_ref(1,:);
            rho0 =  obj.bc(iboun).rho_ref(1,:);
            
            
            j_to_look_at = obj.Iin(ki); %all particles
            % can probably be improved when taking the cell structure into
            % account
            d = distanceToParticles(obj.Xj(j_to_look_at,:),...
                                    obj.Xj(obj.Iin(kb),:),...
                                    outer_normal);
            Imirror_local = d <  obj.kernel_cutoff * max(obj.hj(obj.Iin)); %%
            Imirror       = obj.Iin(ki(Imirror_local));

            %change order (for better readablity , not necessary)
            d = flip(d);
            Imirror_local=flip(Imirror_local);
            Imirror = flip(Imirror);
            

            %%
             m_factor = ones(size(Imirror,1),1);
             % to examine: %tweak on the mass
             tweakmass = obj.exp_settings.tweakmass;
             if tweakmass
                if obj.vj(obj.Iin(end),1)>=0  
                    %linear 
                    mV_factor_min = 0.7;
                    mV_factor_max = 1;
                    dmax = max(d(Imirror_local));
                    drel = (dmax - d(Imirror_local))/dmax;
                    m_factor = (1-drel)* mV_factor_min +...
                                 drel   * mV_factor_max;
                else
                    mV_factor_min = 1.5;
                    mV_factor_max = 1;
                    dmax = max(d(Imirror_local));
                    drel = (dmax - d(Imirror_local))/dmax;
                    m_factor = (1-drel) * mV_factor_min +...
                                 drel    * mV_factor_max;
                end
             end
            % or no adjustment    
            
            if ~isempty(Imirror)
                obj.mirrorParticles(Imirror,d(Imirror_local),outer_normal,rho0,v0,p0,e0,m_factor,v_factor,p_factor); 
            else
                keyboard
            end
        end %mirror particles        
        %%
        function d = point_to_line(~,pt, v1, v2)
              if size(pt,2) == 1 %1d
                  d = abs(pt-v1);
              else
                  Np = size(pt,1);
                  a = [v1 - v2,0];
                  b = [pt - ones(Np,1)*v2,zeros(Np,1)];
                  d = sum(cross(ones(Np,1)*a,b).^2,2).^0.5 ./ norm(a);  
              end
        end
        %%
        function BCnrp (obj, iboun) %copy particles and use characteristic method
            Akb = find(obj.bc(iboun).kb);
            n = obj.bc(iboun).outer_normal;
%             dx = obj.Vj(Akb);
            for i = 1:size(Akb,1)
                kb=Akb(i);
                dx =obj.Vj(kb).^(1/obj.dim);
                rho_ref = obj.bc(iboun).rho_ref(i); 
                v_ref   = obj.bc(iboun).v_ref(i,:);
                p_ref   = obj.bc(iboun).p_ref(i);
                
                %create new particles:
                NN = 5; %amount of particles
                Ig = size(obj.Xj,1)+(1:NN)'; %indice
                copyData(obj, kb, Ig)
                % set coordinate
                x = obj.Xj(kb,:);
                obj.Xj(Ig,:) = ones(NN,1)*x + (1:NN)'*dx*n;
                %---
  
                    
                % possibility 1:
                n = obj.bc(iboun).outer_normal;
                % possibility 2:
%                     n = -obj.bc(iboun).outer_normal;

                J1 = -obj.cj(kb).^2 .*(obj.rhoj(kb) - rho_ref)...
                    + (obj.pj(kb) - p_ref);
                %J1 = 0 always
                J2 = obj.rhoj(kb).*obj.cj(kb) .* (obj.vj(kb,:) - v_ref)*n'...
                    + (obj.pj(kb) - p_ref);

                J3 = - obj.rhoj(kb).*obj.cj(kb) .* (obj.vj(kb,:) - v_ref)*n'...
                + (obj.pj(kb) - p_ref);

%                     if norm(J1)>1e-4
%                         keyboard
%                     end
                % possibility 1:
                J3=0;         
                % possibility 2:
%                     J1=0;
%                     J2=0;


                obj.rhoj(kb) = rho_ref + ...
                    1./(obj.c0j(kb).^2).*(-J1+0.5*J2+0.5*J3);
                obj.vj(Ig,:) = ones(size(Ig,1),1)*(v_ref + ...
                    (1./(2*obj.c0j(kb).*obj.rho0j(kb)) * (J2 - J3))*n);
                obj.pj(Ig,:) = p_ref+...
                    0.5*(J2+J3);                

            end

        end
        %%
        function comp_BCdamping(obj) %todo
           for boun = obj.bc
               if ~isempty(boun.damping_area)
                   error('todo');
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
%                    obj.Fj_tot(kb,:) = (obj.Fj_tot(kb,:).*f_Q + Q.*damp ) ;
% 
% 
%                    %% stop particels
%         %            kb_stop = logical ((obj.Xj>=xb(2)));
%         %            obj.vj(kb_stop) = 0;
               end
           end
          end   
        %%
        function BCnoflow (obj,i)
            boun = obj.bc(i);
            outer_normal = boun.outer_normal;
            v0   = zeros(1,obj.dim);
            e0   = 0;
            p0   = 0;
            rho0 = 0;
            m_factor = 1;
            p_factor  = 1;

            if obj.dim ==1
                v_factor = -1;
            else
               if outer_normal(1)~=0 && outer_normal(2) ==0
                    v_factor = [-1,1];
               elseif outer_normal(1)==0 && outer_normal(2) ~=0
                    v_factor = [1,-1];
               else
                  error('other normal-directions are not implemented right now');
               end
            end

            j_to_look_at = obj.Iin; %all particles
            pXj = obj.Xj(j_to_look_at,:);
            
            
            NN=size(pXj,1);
            bp=boun.bp;
            d = abs(sum((pXj- (ones(NN,1)*bp)) .* (ones(NN,1)*outer_normal),2));
%             if obj.dim == 1
%                  d = abs(pXj-boun.p1);
%             else
%                  d = obj.point_to_line(pXj,boun.p1,boun.p2);
%             end
            Imirror_local = d <  obj.kernel_cutoff * max(obj.hj(obj.Iin));% * 12;
            Imirror       = obj.Iin(Imirror_local);                
            
            if ~isempty(Imirror)
                obj.mirrorParticles(Imirror,d(Imirror_local),outer_normal,rho0,v0,p0,e0,m_factor,v_factor,p_factor);
                obj.comp_pressure();
            end                        
        end
        %%        
        function resetData(obj)
             %reset to initial
            obj.Ighost = [];
            obj.Xj     = obj.Xj(obj.Iin,:);
            obj.vj     = obj.vj(obj.Iin,:);
            
            obj.c0j     = obj.c0j  (obj.Iin,:);
            obj.rho0j   = obj.rho0j(obj.Iin,:);
            if obj.compGhost
                obj.vj_half   = obj.vj_half(obj.Iin,:);
                obj.ej_half   = obj.ej_half(obj.Iin,:);
                obj.rhoj_half = obj.rhoj_half(obj.Iin,:);
            end
            
            obj.ej   = obj.ej(obj.Iin,:);
            obj.mj   = obj.mj(obj.Iin);
            obj.Vj   = obj.Vj(obj.Iin);
            obj.rhoj = obj.rhoj(obj.Iin);
            obj.hj   = obj.hj(obj.Iin);
            obj.cj   = obj.cj(obj.Iin);
            obj.pj   = obj.pj(obj.Iin);
            obj.MG_Gammaj = obj.MG_Gammaj(obj.Iin);
            obj.MG_Sj = obj.MG_Sj(obj.Iin);
            obj.N    = size(obj.Xj,1);
            obj.Imaterial_with_boun = obj.Imaterial;     
            obj.Icomp = obj.Iin;
        end
        %%
        function mirrorParticles (obj,Imirror,d,outer_normal,rho0,v0,p0,e0,m_factor,v_factor,p_factor)            
            %add boundaries:            

            NN = size(Imirror,1);                   
            obj.Ighost = [obj.Ighost;...
                          (obj.N+(1:NN))'];%Imirror;
             
            %copy data
            obj.Xj (obj.N+(1:NN),:) = ...
                      obj.Xj(Imirror,:) + 2*d*outer_normal; 

            obj.vj (obj.N+(1:NN),:) = ...    -...
                      2*ones(NN,1)*v0+obj.vj(Imirror,:).*(ones(NN,1)*v_factor);
                  % ----
            obj.c0j (obj.N+(1:NN),:)   = obj.c0j(Imirror,:);
            obj.rho0j (obj.N+(1:NN),:) = obj.rho0j(Imirror,:);
            if obj.compGhost
                obj.vj_half (obj.N+(1:NN),:) = ...    -...
                          2*ones(NN,1)*v0+obj.vj_half(Imirror,:).*(ones(NN,1)*v_factor);
                obj.ej_half (obj.N+(1:NN),:) = ...    -...
                          2*e0+p_factor.*obj.ej_half(Imirror,:);
                obj.rhoj_half (obj.N+(1:NN),:) = 2*rho0+p_factor.*obj.rhoj_half(Imirror,:);
            end
                   %-----
            obj.ej (obj.N+(1:NN),:) = ... obj.ej(Imirror)
                        2*e0+p_factor.*obj.ej(Imirror);                          
            obj.pj (obj.N+(1:NN),:) = ...
                2*p0+p_factor.*zeros(NN,1).*obj.pj(Imirror);%...                            

            obj.mj (obj.N+(1:NN),:) = m_factor.*obj.mj(Imirror);

            obj.Vj (obj.N+(1:NN),:) = obj.Vj(Imirror);

            obj.hj (obj.N+(1:NN),:) = ...
                      obj.hj(Imirror);
            obj.cj (obj.N+(1:NN),:) = ...
                      obj.cj(Imirror); 
            obj.rhoj (obj.N+(1:NN),:) = ...
                      obj.rhoj(Imirror);
            obj.MG_Gammaj (obj.N+(1:NN),:) = ...
                        obj.MG_Gammaj(Imirror);
            obj.MG_Sj (obj.N+(1:NN),:) = ...
                obj.MG_Sj(Imirror);
            obj.N = size(obj.Xj,1);

            obj.Imaterial_with_boun=[obj.Imaterial_with_boun;
                                     obj.N-NN+1 , obj.N];                                         
                                 
            if obj.compGhost
                obj.Icomp = [obj.Icomp; obj.Ighost];                
            end            
        end 
        %%
        function copyData(obj, Asource, Adestination) %everything is a colum vector
            obj.c0j   (Adestination,:) = obj.c0j(Asource,:);
            obj.rho0j (Adestination,:) = obj.rho0j(Asource,:);
            NN = size(Adestination,1);
            if obj.compGhost
                obj.vj_half (Adestination,:) = ones(NN,1)*obj.vj_half(Asource,:);
                obj.ej_half (Adestination,:) = obj.ej_half(Asource,:);
                obj.rhoj_half (Adestination,:) = obj.rhoj_half(Asource,:);
            end
            obj.Xj (Adestination,:) = ones(NN,1)*obj.Xj(Asource,:);
            obj.vj (Adestination,:) = ones(NN,1)*obj.vj(Asource,:);         
            obj.ej (Adestination,:) = obj.ej(Asource);                          
            obj.pj (Adestination,:) = obj.pj(Asource);                           
            obj.mj (Adestination,:) = obj.mj(Asource);
            obj.Vj (Adestination,:) = obj.Vj(Asource);
            obj.hj (Adestination,:) = obj.hj(Asource);
            obj.cj (Adestination,:) = obj.cj(Asource); 
            obj.rhoj(Adestination,:)= obj.rhoj(Asource);            
            
            NN = size(Adestination,1);
            obj.N = size(obj.Xj,1);
            obj.Imaterial_with_boun=[obj.Imaterial_with_boun;
                                     obj.N-NN+1 , obj.N];                                                                          
            if obj.compGhost
                obj.Icomp = [obj.Icomp; Adestination];                
            end             
        end 
        %%
        function checkDataOnline(obj)          
            %%---------------------
            function checkIfInDomain(obj)
                if (any(any(obj.Xj-  ones(obj.N,1)* obj.Omega(:,1)'<0)) || any(any( ones(obj.N,1)*obj.Omega(:,2)' - obj.Xj <0)))
                    error('some points are out of the domain');
                end
            end  
            %%---------------------
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
            %amount of neibours:
            if ~isempty(obj.AedgesXj)
                nj=sum(abs(obj.AedgesXj(:,obj.active_k_pij)),2);
                disp (['#neigbhours: min=',num2str(min(nj)),...
                                  ', max=',num2str(max(nj)),...
                                  ', mean=',num2str(mean(nj))]);
            end
        end       

    end
end
