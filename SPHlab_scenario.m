% SPH for HVI  - Markus Ganser - TU/e - 2015
classdef SPHlab_scenario < handle
    % This class defines basically the setup of a simulation
    %   - geometry - provides function to generate the geometry
    %   - defines initial state of the particles
    %   - defines plotting adjustments
    % See default values in the properties section.
    % Main functions:
    % - add_geometry
    %   to define the geometry and its material property
    % - add_bc
    %   defines the boundary treatment
    % - create_geometry
    %   creates with add_(.) functions the initial point set and assign the
    %   according values. Furthermoe, it takes care about the BC and if
    %   necessary, it remove particles behind a boundary.
    
    properties

        %% simulation parameter
        dtfactor    % savetyfactor for timestepping with CFL      
        dt          % timestep (if empty, use CFL)
        tend        % Simulation time
        tpause       % time to make a pause the simulation
        eta         % h=eta*dx
        kernel      % M4 | Gauss | Wendland
        kernel_cutoff  % r/h<=kernel_cutoff

        scheme      % m | v
        EOS         % ISO (isothermal) | MG (Mie-Gruneisen) | Water | 
        h_const     % is h constant (true) or dependent on density (false) 
        normalizeOmega % only m-scheme right now
        compOmegaj  % bool - compute correction factor Omega_j (comming from the rho dependency of h)
        %% geometry
        Ntot        % desired amount of particles
        equalmass   % boolean, else equal volume - for initialization (some fluctuation can appear in generation process);
        Omega       % domain  [x_left, x_right]^n
        
        
        geo          % geometry properties of each object (struct)
        Xj           % coordinates
        geo_noise    % some noise in the position
        vj           % velocity                
        Vj           % volume per particle
        ej           % specific internal energy
        
        %% indices domain
        Iin  
        Imaterial  % [1st indice of n-st material, last indice of n-st material]^n
        
        %bc
        damping_area % [lower-left point, upper-right point]^n
        mirrorParticlesj
        bc;       %'p1'; point on the boundary line 
                  %'outer_normal'
                  %'damping_area' (experimantal)
        
        %% material parameter
        rho0j      % density
        rhoj
        c0j       % speed of sound
        cj
        MG_Gammaj  % MieGruneisen parameter
        MG_Sj      % relation between shock and particle velocity
        beta      % for surface tension (scalar)
        mu        % for dissipation
        art_diss_para  % parameters for artificial dissipation (alpha_mass, alpha_viscosity, beta_viscosity, alpha viscosity)
        
        
        % experimantal settings (struct)
        exp_settings
        
        % 
        g_ext    %gravity

        %exact solution
        exact_sol
        %% IO
        % output
        write_data
        save_dt             %timestep for saving the data
        output_name        
        save_as_figure %bool
        figure_format
        % some plotting properties
        
        plot_dt         % plotting timestep
        plot_style       % 1D - scatter  
                         % 2D - scalar: trisurf | patches; field: quiver
        plot_quantity    %which quantity shall be plottet: eg 'xe' for position and energy
         % v-velocity, x-position, p-pressure, d-density, f-forces, m-massflux, e-energy
        fixaxes         %struct to define the axes
        
        plotconfig      % struct: see members below
        %%        
        Neval           % amount of evaluation points 
        Nss             % amount of supersampling points (in each direction)
        %movie settings                        
        save_as_movie
        movie_name
        
        iter   % iterator - what indice would a new particle have        
        iter_boun %counts the boundary conditions
        iter_geo  %counts the different geometries
    end
    
    methods
        % Constructor
        function obj = SPHlab_scenario(filename)
           % some standard parameter 
           obj.kernel = 'M4'; 
           obj.kernel_cutoff = 2;
           obj.scheme = 'm';
           obj.EOS    = 'ISO';
           obj.equalmass = false;
           obj.h_const   = false;
           obj.normalizeOmega = false;
           obj.compOmegaj= false;
           obj.dt        = [];
           obj.dtfactor  = 0.2;
           obj.tpause    = inf;

           obj.beta = 0; %material dependent!
           obj.mu   = 0;
           obj.art_diss_para = struct('alpha_mass',0.5,... %Iason2014 (Spheric)
                                      'beta_mass',0,...
                                      'alpha_viscosity',1,...
                                      'beta_viscosity',2,...
                                      'alpha_energy',0.3,...
                                      'beta_energy',0);
           obj.geo = struct([]);

           obj.geo_noise = 0;
           obj.g_ext = []; 

           %IO
           obj.plot_dt = inf;
           obj.plot_quantity = 'xp';
           obj.plot_style = struct('x','scatter',...              
                                   'p','patches',...
                                   'd','patches',...
                                   'm','patches',...
                                   'v','quiver',...
                                   'f','quiver',...
                                   'e','patches',...
                                   'c','patches');
           obj.plotconfig = struct();
           obj.plotconfig.drawfluxes = true;            
           obj.plotconfig.latexplot  = false;
           obj.plotconfig.transpose  = false; %transpose the subplots
           obj.plotconfig.figuresize = []; %empty: make default
           obj.plotconfig.figurename = [];
           
           obj.Neval = 0;
           obj.Nss   = 1; %(no supersampling)
           
           obj.save_dt = 0;
           obj.save_as_movie = false;
           obj.save_as_figure = false;
           obj.figure_format = 'eps';
           obj.movie_name = 'out';
           obj.write_data = false; 
           obj.output_name ='data/data_out';
           obj.fixaxes = struct('x',[],'v',[],'p',[],'d',[],'m',[],...
                                'f',[],'e',[],'c',[]);
           obj.bc      = struct([]);
           obj.exp_settings = struct(...
                                     'tweakmass',false...  
                                     );
                                 
           obj.exact_sol =[];
           %some initialization
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.damping_area = [];
           obj.mirrorParticlesj = [];
           obj.iter  = 1; 
           obj.iter_boun = 1;
           obj.iter_geo = 1;

           if nargin  > 0
               %% read file and overwrite some properties
               % with excisting scenario:
               filename_scen = [filename,'_scen.mat'];
               if exist (filename_scen,'file')
                   load(filename_scen);
                   obj = obj_scen;  %overwrite everything  
               else
                   %without (e.g. from LimeSPH)
                   group = '/0';
                   obj_IO = SPHlab_IO();                   
                   obj_IO.read_hdf5(obj,filename,group);
                   warning('Make sure that all properties (like Omega) are defined!');
               end
               obj.write_data = false; 
               obj.output_name =[obj.output_name,'_new'];
           end
           
           
           disp('---------------------------------');
           disp('-------------SPHlab--------------');
           disp('---------------------------------');
        end   
        
        %%
        function name = get_movie_name(obj)
            movie_dir='movies/'; %ToDo: anders machen
            movie_format='.mp4';
            name=[movie_dir,obj.movie_name,movie_format];
        end
        %%  
        function dispdata(obj)
           disp(obj)             
        end
        %%
        function set_kernel(obj,kernel)
           obj.kernel = kernel;
           if strcmp(kernel, 'Wendland')
               obj.kernel_cutoff = 2;
           elseif strcmp(kernel, 'M3') %after violeau (M4 in Price)
               obj.kernel_cutoff = 2;
           elseif strcmp(kernel, 'M4')
               obj.kernel_cutoff = 2.5;
           elseif strcmp(kernel, 'Gauss')
               obj.kernel_cutoff = 4;
           elseif strcmp(kernel, 'linear')
               warning('Linear kernel is experimental!');
               obj.kernel_cutoff = 2;
           else
               error('Kernel not supported!');
           end
        end
        %%
        function add_geometry(obj,omega_geo, rho0, v0, c0,e0,MG_Gamma,MG_S,Nfactor)
            if nargin <6
                e0 = 0;
            end     
            if nargin <7                
                if strcmp(obj.EOS,'MG')
                    error('Please define some Mie-Gruneisen parameters!');
                else
                    MG_Gamma = 0; %just set something
                    MG_S = 0;
                end
            end
            
            if nargin <9
                Nfactor =1;
            end
            obj.geo(obj.iter_geo).omega_geo = omega_geo;
            obj.geo(obj.iter_geo).rho0      = rho0;
            obj.geo(obj.iter_geo).v0        = v0;
            obj.geo(obj.iter_geo).c0        = c0;
            obj.geo(obj.iter_geo).e0        = e0;
            obj.geo(obj.iter_geo).MG_Gamma  = MG_Gamma;
            obj.geo(obj.iter_geo).MG_S      = MG_S;

            obj.geo(obj.iter_geo).Nfactor   = Nfactor;
            
            obj.iter_geo = obj.iter_geo +1;
        end
        %%
        function create_geometry(obj)
            %if only a total amount of particle is given, try to distribute
            %the particles in a fair way
            dim = size(obj.Omega,1);
            if length(obj.Ntot) == 1 
                m_tot = 0;
                V_tot = 0;
                %compute overall volume and mass
                for i = 1:size(obj.geo,2)
                   if dim == 2 && size(obj.geo(i).omega_geo,1) == 1 %circle
                     obj.geo(i).V = pi*obj.geo(i).omega_geo(1).^2;
                   else %usual rectangular
                     obj.geo(i).V = prod(diff(obj.geo(i).omega_geo'));
                   end
                   V_tot = V_tot + obj.geo(i).V;
                   obj.geo(i).m = obj.geo(i).V * obj.geo(i).rho0;
                   m_tot = m_tot + obj.geo(i).m;
                end
                %compute particels per geometry object
                for i = 1:size(obj.geo,2)
                    if obj.equalmass
                        obj.geo(i).N = obj.Ntot/m_tot * obj.geo(i).m;
                    else
                        obj.geo(i).N = obj.Ntot/V_tot * obj.geo(i).V;
                    end
                    %adjust amount of particles
                    obj.geo(i).N = round(obj.geo(i).Nfactor * obj.geo(i).N);
                end
            % otherwise use predefined number
            else
                if length(obj.Ntot) ~= obj.iter_geo-1
                   error ('Size of Ntot does not fit to the number of geometries.');                     
                end
                for i = 1:size(obj.geo,2)
                   obj.geo(i).N = obj.Ntot(i);                     
                end
            end                        
            
            for i = 1:size(obj.geo,2)                
                %create geometry
                if dim == 1
                    [x,Vparticle] = add_line(obj,obj.geo(i).omega_geo, obj.geo(i).N);                    
                elseif dim == 2
                    if size(obj.geo(i).omega_geo,1) == 1 % circle
                        [x,Vparticle] = add_circle(obj,obj.geo(i).omega_geo, obj.geo(i).N);
                    else %rectangle                        
                        hexa=false; %make hexagonal particle distribution
                        [x,Vparticle] = add_rectangle(obj,obj.geo(i).omega_geo, obj.geo(i).N,hexa);
                    end
                else
                    error('only dim=1,2 are supported');
                end
                
                % cut on boundary, remove particles:
                for boun = obj.bc
                   if strcmp(boun.type,'nrd')
                       continue,
                   end
                   NN=size(x,1);
                   bp=boun.bp;
                   n =boun.outer_normal;
                   %compute distance to boundary point (in normal
                   %direction)
                   d = sum((x- (ones(NN,1)*bp)) .* (ones(NN,1)*n),2);
                   % negativ distance is ok, delete particles with positive
                   % distance
                   Igood = d<0;
                   x=x(Igood,:);
                   if any(d>=0)
                      disp('cut on boundary - remove particles') 
                   end                  
                end                
                
                %update indice iterator
                Np= size(x,1);
                I = (obj.iter : obj.iter+Np-1)';
                obj.iter = obj.iter+Np;
                
                %set properties:
                obj.Xj(I,:)    = x;
                obj.Vj(I,:)    = Vparticle;
                obj.vj(I,:)    = ones(size(I,1),1)*obj.geo(i).v0;     
                obj.c0j(I,1)   = obj.geo(i).c0;
                obj.cj(I,1)    = obj.geo(i).c0;
                obj.ej(I,1)    = obj.geo(i).e0;
                obj.MG_Gammaj(I,1)    = obj.geo(i).MG_Gamma;
                obj.MG_Sj(I,1)    = obj.geo(i).MG_S;
                obj.rho0j(I,1) = obj.geo(i).rho0;
                obj.rhoj(I,1)  = obj.geo(i).rho0;
                obj.Imaterial = [obj.Imaterial;...
                                 [I(1) I(end)] ];
                obj.Iin   = [obj.Iin; I];
                disp(['Object ',num2str(i),' contains ',num2str(size(I,1)), ' particles.']);
            end
        end
        %%
        function add_bc(obj,type,bp,outer_normal,damping_area) %type: nr-c,nr-m,nr-p,noflow
            if nargin < 5
                damping_area=[];
            end
            obj.bc(obj.iter_boun).type = type;
            obj.bc(obj.iter_boun).outer_normal = outer_normal./norm(outer_normal);
            obj.bc(obj.iter_boun).bp = bp; %boundary point 

            obj.bc(obj.iter_boun).damping_area = damping_area;
            obj.iter_boun = obj.iter_boun + 1;           
        end
        
        %% % some geometry functions % %%
        
        %% generic function to create a rectangel in 1d and 2d
        function [xy,dxy] = rectangle (~,omega_geo, N,hexa)
            if nargin<4
                hexa=false;
            end
           %generate N points in omega_geo
            len_xy  = diff(omega_geo');
            dim     = size(omega_geo,1);
            if (dim == 0)
                error('something wrong with the dimension')
            end
            
            
            len_particle   = (prod(len_xy)/N)^(1/dim);
            % amount of particle needed for an equidistant mesh
            N_xy = len_xy/len_particle; 
            % round this number            
            N_xy = round(N_xy);
            Np = prod(N_xy);
            dxy  = len_xy ./ N_xy;    %mesh size    
            % ||dx/2| X |dx| X |dx| X ... X |dx| X |dx/2|| 
            
            for i = 1:dim
                xtemp = (linspace(omega_geo(i,1)+dxy(i)/2,...
                                  omega_geo(i,2)-dxy(i)/2,...
                                  N_xy(i)))';                
                if i == 2 %add second dimension
                    % create mesh
                    [t1,t2]= ndgrid(xy,xtemp);
                     
                    if hexa
                        %hexa:
%                         %move all a bit to the left
%                         t1=t1-0.25*dxy(1);
%                         %move every second row to the right
%                         I=1:2:size(t1,2);
%                         t1(:,I)=t1(:,I)+0.5*dxy(1);
%                         
                        
                        %hexa: (2nd version vertical shift)
                        %move all a bit to the bottom
                        t2=t2-0.25*dxy(2);
                        %move every second row to the top
                        I=1:2:size(t2,1);
                        t2(I,:)=t2(I,:)+0.5*dxy(2);                        
                    end
                    xx = reshape(t1,[1,Np]); 
                    yy = reshape(t2,[1,Np]);
                    xy=[xx',yy'];                    
                else
                    xy = xtemp;
                end
            end
        end        
        %%
        function [xy,Vparticle] = add_rectangle(obj,omega_geo, N,hexa) %for 1 and 2 dimension
            if nargin<4
                hexa=false;
            end
            [xy,dxy] = rectangle (obj,omega_geo, N,hexa);
            % Volume of one particle
            Vparticle = prod(dxy);          
        end
        %%
        function [xy,Vparticle] = add_circle(obj,omega_geo, N) %create circle in an easy way 
            %first, create an rectangle
            r=omega_geo(1);
            c=omega_geo(2:3);
            omega_rec = [c(1)-r,c(1)+r;c(2)-r,c(2)+r];
            [xy,dxy] = rectangle (obj,omega_rec, N);
            % Volume of one particle
            Vparticle = prod(dxy);          
            %remove particles outside of this circle
            I=(sum((xy-(ones(size(xy,1),1)*c)).^2,2).^0.5)<r;
            xy=xy(I,:);
        end
        %%           
        function [xy,Vparticle] = add_line(obj,omega_geo, N)
            %use more general function
            [xy,Vparticle] = add_rectangle(obj,omega_geo, N);
        end        
        %%
        function dx = dx_min (obj)
            dim = size(obj.Xj,2);
            dx = min(obj.Vj)^(1/dim);
        end
        function dx = dx_max (obj)
            dim = size(obj.Xj,2);
            dx = max(obj.Vj)^(1/dim);
        end        
                              
    end    
end

