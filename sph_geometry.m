classdef sph_geometry < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % GEOMETRY class - creates the initial geometry 
    
    properties
        % particle properties
        Xj           % coordinates
        V0particle   % Volume of the particle
        vj           % velocity
        iter         % iterator - how much particle are already in Xj
        
    end
    
    methods
        % Constructor
        function obj = sph_geometry()
           obj.Xj    = []; 
           obj.vj    = []; 
           obj.V0particle = [];
           obj.iter  = 1;
        end
        
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

