classdef sph_scenario < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % scenario class - with geometry functions
    
    properties

        %% simulation parameter
        dim         % 1 or 2
        dx          % initial particle distance
        omega       % savetyfactor for timestepping
        tend        % Simulation time
        eta         % h=eta*dx
        eta2        % eta2*h is the cutoff radius 
        kernel      % M4 | gauss
        %% geometry
        Omega       % domain  [x_left, x_right]^n
        obj_geo     % geometry class
        
        Xj           % coordinates
        V0particle   % Volume of the particle
        vj           % velocity                
        
        %% indices domain
        Iboun
        Iin  
        Imaterial  % [1st indice of n-st material, last indice of n-st material]^n
       
        %% material parameter
        rho0     % relative density
        c0       % speed of sound
        beta     % for surface tension
        mu       % for dissipation
        
        % 
        g_ext    %gravity

        %% IO
        % input
        read_file
        input_name
        % output
        
        
        % some plotting properties
        plot_dt         % plotting timestep
        plotstyle       % 1D: x-position, v-velocity, p-pressure, d-density, f-forces
                        % 2D: scatter | trisurf | patches  
        %movie settings                        
        save_as_movie
        movie_name
        
        iter   % iterator - what indice would a new particle have        
    end
    
    methods
        % Constructor
        function obj = sph_scenario()
           obj.Xj    = []; 
           obj.vj    = []; 
           obj.V0particle = [];
           obj.Iboun = [];
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.g_ext = [0,0];
           obj.save_as_movie = false;
           obj.movie_name = 'out';
           obj.plotstyle = 'scatter';
           obj.read_file = false;
           obj.kernel = 'M4'; %standard kernel
           obj.omega  = 0.5;
           obj.beta = 0;
           obj.mu   = 0;
           
           obj.iter  = 1;
        end       
        
        function checkIfAlreadySet(obj)
            if ~isempty(obj.Iin)
                error('A scenario is already set!');
            end
        end
        
        function name = get_movie_name(obj)
            movie_dir='movies/'; %ToDo: anders machen
            movie_format='.avi';
            name=[movie_dir,obj.movie_name,movie_format];
        end
          
        function dispdata(obj)
           disp(obj)             
        end
        
        
        %% % some geometry functions % %%
         function I=add_line1d(obj,Astartpoint, Aendpoint,v0,dx)
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
                    obj.vj = [obj.vj;v0];     %add velocity
                    obj.V0particle = [obj.V0particle; dx];
                    obj.iter=obj.iter+1;
                    x_par = x_par + e_edge*dx;
                end
            end
            second_ind=obj.iter-1;
            I=(first_ind:second_ind)';
        end
        
        function I=add_line2d(obj,Astartpoint, Aendpoint,v0,dx,layer,noisefactor)  
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
                        x_par_temp = x_par + (shift-1)*e_edge*[0,-1;0,1]*dx;
                        noise = dx*noisefactor*(2*rand(1,2)-1); %add some noise
                        obj.Xj = [obj.Xj;x_par_temp+noise];  %add point
                        obj.vj = [obj.vj;v0];     %add velocity
                        obj.V0particle = [obj.V0particle; dx*dx];
                        obj.iter=obj.iter+1;
                    end
                    x_par = x_par + e_edge*dx;
                end
            end
            second_ind=obj.iter-1;
            I=(first_ind:second_ind)';
            
        end
        
        function I=add_rectangle2d(obj,Alowerleftcorner, Aupperrightcorner,...
                v0, dx,dy,noisefactor)
            first_ind=obj.iter;
            for k=1:size(Alowerleftcorner,1)
                lowerleftcorner = Alowerleftcorner(k,:);
                upperrightcorner= Aupperrightcorner(k,:);
                x_par=lowerleftcorner;
                while (x_par(2) <= upperrightcorner(2));
                    while (x_par(1)<= upperrightcorner(1));
                        noise = dx*noisefactor*(2*rand(1,2)-1);
                        obj.Xj = [obj.Xj;x_par+noise];  %add point
                        obj.vj = [obj.vj;v0];     %add velocity
                        obj.V0particle = [obj.V0particle; dx*dy];
                        obj.iter=obj.iter+1;
                        x_par(1)= x_par(1)+dx;
                    end
                    x_par(1)= lowerleftcorner(1);
                    x_par(2)= x_par(2)+dy;
                end  
            end
            second_ind=obj.iter-1;
            I=(first_ind:second_ind)';
        end
        
               
    end    
end

