classdef sph_scenario < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % scenario class - with geometry functions
    
    properties

        %% simulation parameter
        dtfactor    % savetyfactor for timestepping with CFL      
        dt          % timestep (if empty, use CFL)
        tend        % Simulation time
        eta         % h=eta*dx
        kernel      % M4 | Gauss | Wendland
        kernel_cutoff  % r/h<=kernel_cutoff

        scheme      % m | v
        h_const     % is h constant (true) or dependent on density (false) 
        %% geometry
        Ntot        % desired amount of particles
        equalmass   % boolean, else equal volume - for initialization (some fluctuation can appear in generation process);
        Omega       % domain  [x_left, x_right]^n
        
        
        geo          % geometry properties of each object (struct)
        Xj           % coordinates
        geo_noise    % some noise in the position
        vj           % velocity                
        Vj           % volume per particle
        
        %% indices domain
        Iin  
        Imaterial  % [1st indice of n-st material, last indice of n-st material]^n
        
        %bc
        damping_area % [lower-left point, upper-right point]^n
        mirrorParticlesj
        bc;       %'p1';'p2;   for noflow, two points on boudary
                  %'mirrorParticlesj' for no-reflecting
                  %'outer_normal'
                  %'damping_area'
        
        %% material parameter
        rho0j      % density
        rhoj
        c0j       % speed of sound
        cj
        beta      % for surface tension (scalar)
        mu        % for dissipation
        
        % 
        g_ext    %gravity

        %% IO
        % input
        read_data
        input_name
        % output
        write_data
        output_name
       
        % some plotting properties
        plot_dt         % plotting timestep
        plot_style       % 1D;x: scatter  
                         % 2D - scalar: trisurf | patches; field: quiver
        plot_quantity    %which quantity shall be plottet v-velocity, x-position, p-pressure, d-density, f-forces
        fixaxes         %struct to define the axes           
        %movie settings                        
        save_as_movie
        movie_name
        
        iter   % iterator - what indice would a new particle have        
        iter_boun %counts the boundary conditions
        iter_geo  %counts the different geometries
    end
    
    methods
        % Constructor
        function obj = sph_scenario()
           % some standard parameter 
           
           obj.kernel = 'Wendland'; 
           obj.kernel_cutoff = 2;
           obj.scheme = 'm';
           obj.equalmass = false;
           obj.h_const   = false;           
           obj.dt        = [];
           obj.dtfactor  = 0.5;
           
           obj.beta = 0; %material dependent!
           obj.mu   = 0;
                   
           obj.geo = struct([]);
           
           obj.geo_noise = 0;
           obj.g_ext = []; 
           
           %IO
           obj.save_as_movie = false;
           obj.movie_name = 'out';
           obj.plot_quantity = 'xp';
           obj.plot_style = struct('x','scatter',...              
                                   'p','patches',...
                                   'd','patches',...
                                   'v','quiver',...
                                   'f','quiver',...
                                   'e','patches');
           obj.read_data = false;           
           obj.write_data = false; 
           obj.output_name ='data_out.h5';
           obj.fixaxes = struct('x',[],'v',[],'p',[],'d',[],'f',[],'e',[]);
           obj.bc      = struct([]);
           
           %some initialization
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.damping_area = [];
           obj.mirrorParticlesj = [];
           obj.iter  = 1; 
           obj.iter_boun = 1;
           obj.iter_geo = 1;
        end   
        
        %%
        function name = get_movie_name(obj)
            movie_dir='movies/'; %ToDo: anders machen
            movie_format='.avi';
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
           elseif strcmp(kernel, 'M4') %after violeau (M4 in Price)
               obj.kernel_cutoff = 2.5;
           elseif strcmp(kernel, 'Gauss')
               obj.kernel_cutoff = 2;
           else
               error('Kernel not supported!');
           end
        end
        %%
        function add_geometry(obj,omega_geo, rho0, v0, c0,Nfactor)
            if nargin <6
                Nfactor =1;
            end
            obj.geo(obj.iter_geo).omega_geo = omega_geo;
            obj.geo(obj.iter_geo).rho0      = rho0;
            obj.geo(obj.iter_geo).v0        = v0;
            obj.geo(obj.iter_geo).c0        = c0;
            obj.geo(obj.iter_geo).Nfactor   = Nfactor;
            
            obj.iter_geo = obj.iter_geo +1;
        end
        %%
        function create_geometry(obj)
            m_tot = 0;
            V_tot = 0;
            dim = size(obj.Omega,1);
            %compute overall volume and mass
            for i = 1:size(obj.geo,2)
               obj.geo(i).V = prod(diff(obj.geo(i).omega_geo'));
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
            for i = 1:size(obj.geo,2)
                %create geometry
                if dim == 1
                    [x,Vparticle,I] = add_line(obj,obj.geo(i).omega_geo, obj.geo(i).N);                    
                elseif dim == 2
                    [x,Vparticle,I] = add_rectangle(obj,obj.geo(i).omega_geo, obj.geo(i).N);                                        
                else
                    error('only dim=1,2 are supported');
                end
                %set properties:
                obj.Xj(I,:)    = x;
                obj.Vj(I,:)    = Vparticle;
                obj.vj(I,:)    = ones(size(I,1),1)*obj.geo(i).v0;     
                obj.c0j(I,1)   = obj.geo(i).c0;
                obj.cj(I,1)    = obj.geo(i).c0;
                obj.rho0j(I,1) = obj.geo(i).rho0;
                obj.rhoj(I,1)  = obj.geo(i).rho0;
                obj.Imaterial = [obj.Imaterial;...
                                 [I(1) I(end)] ];
                obj.Iin   = [obj.Iin; I];
                disp(['Object ',num2str(i),' contains ',num2str(size(I,1)), ' particles.']);
            end
        end
        %%
        function add_bc_nr(obj,p1,p2,outer_normal) %no reflecting bc
            obj.bc(obj.iter_boun).type = 'nr';
            obj.bc(obj.iter_boun).outer_normal = outer_normal;
            obj.bc(obj.iter_boun).p1 = p1; 
            obj.bc(obj.iter_boun).p2 = p2;
            obj.bc(obj.iter_boun).damping_area = [];
            obj.iter_boun = obj.iter_boun +1;
        end
        %%
        function add_bc_noflow(obj,p1,p2,outer_normal)
            obj.bc(obj.iter_boun).type = 'noflow';
            obj.bc(obj.iter_boun).p1 = p1;
            obj.bc(obj.iter_boun).p2 = p2;
            obj.bc(obj.iter_boun).outer_normal = outer_normal;
            obj.bc(obj.iter_boun).damping_area = [];
            obj.iter_boun = obj.iter_boun +1;            
        end
        
        %% % some geometry functions % %%
        
        %% generic function to create a rectangel in 1 and 2d
        function [xy,dxy] = rectangle (~,omega_geo, N)
           %generate N points in omega_geo
            len_xy  = diff(omega_geo');
            dim     = size(omega_geo,1);
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
                    xx = reshape(t1,[1,Np]);
                    yy = reshape(t2,[1,Np]);
                    xy=[xx',yy'];
                else
                    xy = xtemp;
                end
            end
        end        
        %%
        function [xy,Vparticle,I] = add_rectangle(obj,omega_geo, N) %for 1 and 2 dimension
            [xy,dxy] = rectangle (obj,omega_geo, N);
            % Volume of one particle
            Vparticle = prod(dxy);          
            Np = size(xy,1);
            %update iterator
            I = (obj.iter : obj.iter+Np-1)';
            obj.iter = obj.iter+Np;
        end
        %%           
        function [xy,Vparticle,I] = add_line(obj,omega_geo, N)
            %use more general function
            [xy,Vparticle,I] = add_rectangle(obj,omega_geo, N);
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
        
        
        %% outdated
        function I=add_line1d_old(obj,Astartpoint, Aendpoint)
            first_ind=obj.iter;
            for k=1:size(Astartpoint,1)
                startpoint=Astartpoint(k,:);
                endpoint=Aendpoint(k,:);
                e_edge=endpoint-startpoint;
                l_edge = norm(endpoint-startpoint);
                e_edge = e_edge/l_edge;
                x_par=startpoint;
                while (norm(x_par-startpoint)<=l_edge)
                    obj.Xj = [obj.Xj;x_par];  %add point                 
                    obj.iter=obj.iter+1;
                    x_par = x_par + e_edge*obj.dx;
                end
            end
            second_ind = obj.iter-1;
            I = (first_ind:second_ind)';
        end        
        %% outdated
        function I=add_line2d(obj,Astartpoint, Aendpoint,layer)  
            %place particles from startpoint to endpoint with distance dx
            %and with #layer
            first_ind=obj.iter;
            for k=1:size(Astartpoint,1)
                startpoint=Astartpoint(k,:);
                endpoint=Aendpoint(k,:);
                e_edge=endpoint-startpoint;
                l_edge = norm(endpoint-startpoint);
                e_edge = e_edge/l_edge;
                x_par=startpoint;
                while (norm(x_par-startpoint)<=l_edge)
                    for shift = 1:layer;
                        x_par_temp = x_par + (shift-1)*e_edge*[0,-1;0,1]*obj.dx;
                        noise = obj.dx*obj.geo_noise*(2*rand(1,2)-1); %add some noise
                        obj.Xj = [obj.Xj;x_par_temp+noise];  %add point
                        obj.iter=obj.iter+1;
                    end
                    x_par = x_par + e_edge*obj.dx;
                end
            end
            second_ind=obj.iter-1;
            I=(first_ind:second_ind)';            
        end
        %% outdated
        function I=add_rectangle2d(obj,Alowerleftcorner, Aupperrightcorner)
            first_ind=obj.iter;
            for k=1:size(Alowerleftcorner,1)
                lowerleftcorner = Alowerleftcorner(k,:);
                upperrightcorner= Aupperrightcorner(k,:);
                x_par = lowerleftcorner;
                while (x_par(2) <= upperrightcorner(2));
                    while (x_par(1)<= upperrightcorner(1));
                        noise = obj.dx*obj.geo_noise*(2*rand(1,2)-1);
                        obj.Xj = [obj.Xj;x_par+noise];  %add point
                        obj.iter=obj.iter+1;
                        x_par(1)= x_par(1)+obj.dx;
                    end
                    x_par(1)= lowerleftcorner(1);
                    x_par(2)= x_par(2)+obj.dx;
                end  
            end
            second_ind=obj.iter-1;
            I=(first_ind:second_ind)';
        end
        %% 
        function addproperties(obj, I, Vp, rho0, v0,c0) %constant mass/Volume
            m0 = Vp * rho0;  
            obj.mj(I,1)    = m0;
            obj.vj(I,:)    = ones(size(I,1),1)*v0;     
            obj.c0j(I,1)   = c0;
            obj.cj(I,1)    = c0;
            obj.rho0j(I,1) = rho0;
            obj.rhoj(I,1)  = rho0;
            obj.Imaterial = [obj.Imaterial;...
                             [I(1) I(end)] ];
            obj.Iin   = [obj.Iin; I];

        end       
        %%
    end    
end

