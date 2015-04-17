classdef sph_scenario < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % scenario class - with geometry functions
    
    properties

        %% simulation parameter
        dx          % initial particle distance
        dtfactor    % savetyfactor for timestepping
        tend        % Simulation time
        eta         % h=eta*dx
        eta2        % eta2*h is the cutoff radius 
        kernel      % M4 | Gauss | Wendland
        scheme      % m | v
        h_const     % is h constant (true) or dependent on density (false) 
        %% geometry
        Omega       % domain  [x_left, x_right]^n
        
        Xj           % coordinates
        geo_noise    % some noise in the position
        V0particle   % Volume of the particle
        vj           % velocity                
        mj           % mass
        
        %% indices domain
        Iboun
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
        beta      % for surface tension
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
        plotstyle       % 1D: x-position, v-velocity, p-pressure, d-density, f-forces
                        % 2D: scatter | trisurf | patches  
        fixaxes         %struct to define the axes           
        %movie settings                        
        save_as_movie
        movie_name
        
        iter   % iterator - what indice would a new particle have        
        iter_boun %counts the boundary conditions
    end
    
    methods
        % Constructor
        function obj = sph_scenario()
           % some standard parameter 
           obj.kernel = 'Wendland'; 
           obj.scheme = 'm';
           obj.h_const = false;           
           obj.dtfactor  = 0.5;
           obj.beta = 0;
           obj.mu   = 0;
           obj.geo_noise = 0;
           obj.g_ext = [0,0]; % not in use yet
           
           %IO
           obj.save_as_movie = false;
           obj.movie_name = 'out';
           obj.plotstyle = 'scatter';
           obj.read_data = false;           
           obj.write_data = false; 
           obj.output_name ='data_out.h5';
           obj.fixaxes = struct('x',[],'v',[],'p',[],'d',[],'f',[]);
           obj.bc      = struct([]);
           %some initialization
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.damping_area = [];
           obj.mirrorParticlesj = [];
           obj.iter  = 1; 
           obj.iter_boun = 1;
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
        function addproperties(obj, I, Vp, rho0, v0,c0) %constant mass/Volume
            m0 = Vp * rho0;            
            obj.vj(I,:)    = ones(size(I,1),1)*v0;     
            obj.c0j(I,1)   = c0;
            obj.cj(I,1)    = c0;
            obj.rho0j(I,1) = rho0;
            obj.rhoj(I,1)  = rho0;
            obj.mj(I,1)    = m0;
            obj.Imaterial = [obj.Imaterial;...
                             [I(1) I(end)] ];
            obj.Iin   = [obj.Iin; I];

        end
        %%
        function add_bc_nr(obj,mirrorParticlesj,outer_normal)
            obj.bc(obj.iter_boun).mirrorParticlesj = mirrorParticlesj;
            obj.bc(obj.iter_boun).outer_normal = outer_normal;
            obj.bc(obj.iter_boun).p1 = []; %just bookkeeping
            obj.bc(obj.iter_boun).p2 = [];
            obj.bc(obj.iter_boun).damping_area = [];
            obj.iter_boun = obj.iter_boun +1;
        end
        %%
        function add_bc_noflow(obj,p1,p2,outer_normal)
            obj.bc(obj.iter_boun).p1 = p1;
            obj.bc(obj.iter_boun).p2 = p2;
            obj.bc(obj.iter_boun).outer_normal = outer_normal;
            obj.bc(obj.iter_boun).mirrorParticlesj = [];
            obj.bc(obj.iter_boun).damping_area = [];
            obj.iter_boun = obj.iter_boun +1;            
        end
        
        %% % some geometry functions % %%
        %%
        function I=add_line1d(obj,Astartpoint, Aendpoint)
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
        %% 
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
        %%
        function I=add_rectangle2d(obj,Alowerleftcorner, Aupperrightcorner)
            first_ind=obj.iter;
            for k=1:size(Alowerleftcorner,1)
                lowerleftcorner = Alowerleftcorner(k,:);
                upperrightcorner= Aupperrightcorner(k,:);
                x_par=lowerleftcorner;
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
        
        
               
    end    
end

