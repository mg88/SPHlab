% SPH for HVI  - Markus Ganser - TU/e - 2015
% main function

function sph
    %startup
    close all; 
    obj_scen = sph_scenarios();
    %% set test scenario
    %set_1d_boun(obj_scen);
    set_1d_riemann(obj_scen);
    %set_2d_riemann(obj_scen);
    %set_2d_moving_lines(obj_scen);
    %set_2d_impact(obj_scen);
    %set_2d_two_squares(obj_scen);
    %set_2d_square(obj_scen);
    
    %% create particle class
    obj_particles=sph_particles(obj_scen);

    %% movie
    mfigure=figure;
    if obj_scen.save_as_movie
        set(gca,'DataAspectRatio',[1,1,1]);
        movie_name = get_movie_name(obj_scen);
        vidObj     = VideoWriter(movie_name);
        open(vidObj);
    end 
 %% time iteration
    ttic = tic;
    t    = 0;
    dt = obj_scen.dt;
    while t < obj_scen.tend
       % tic
        update_full_step(obj_particles,dt)
        comp_pressure(obj_particles, obj_scen.Ca);
        comp_volume(obj_particles);
        search_neighbours(obj_particles);
        update_kernel(obj_particles)
        comp_forces(obj_particles,obj_scen.mu,obj_scen.beta)
        comp_dRho(obj_particles)
        update_half_step(obj_particles,dt)
        update_postion(obj_particles,dt)
        %% plotting
        if (mod(t,obj_scen.plot_dt)<dt)
            plot_data(obj_particles,t,obj_scen.plotstyle,mfigure);
            if obj_scen.save_as_movie
                currFrame = getframe(mfigure);
                writeVideo(vidObj,currFrame);
            end
        end
        t=t+dt;
        %toc
    end
    toc(ttic)
    if obj_scen.save_as_movie
      close(vidObj);
      disp (['movie saved as ',movie_name]);
   end
end

