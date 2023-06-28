for xy = 1:5

clear all
close all

%%%%%%%%%%%%%%%%%%%%%% PARTICLES NUMBER %%%%%%%%%%%%%%%%%%%%%%%%
num_part        = 1;%number of particles after division

%%%%%%%%%%%%%%%%%%%%%% FILLING TABLES %%%%%%%%%%%%%%%%%%%%%%%%%%
posx_part1      = zeros(num_part,1);posx_part2 = zeros(num_part,1);
posy_part1      = zeros(num_part,1);posy_part2 = zeros(num_part,1);
posx_circle     = zeros(num_part,1);posy_circle = zeros(num_part,1);
posy_virtual    = zeros(num_part,1);posx_virtual = zeros(num_part,1);
vect_force_virt_x = zeros(num_part,1);vect_force_virt_y= zeros(num_part,1);
amp_force_virt  = zeros(num_part,1);

dist_part       = zeros(num_part,num_part);
vect_force_x    = zeros(num_part,num_part);
vect_force_y    = zeros(num_part,num_part);
amp_force       = zeros(num_part,num_part);
dist_virtual    = zeros(num_part,num_part);

%%%%%%%%%%% GEOMETRY AND INITIAL POSITIONS OF PARTICLES %%%%%%%%
domain_radius = 40;%domain of the circle

%index_random if you want to have the particle taken at random places in
%the circle
index_random = 0;
%%geometric placement of initial positions, you can do random or use fixed
%%sets

%initial sets of positions for two particles if we choose their position
%never start excatly on 0
set = 3;
if set == 1
     px1 = 0;
     py1 = 0.98*domain_radius;
     px2 = 0.3*domain_radius;
     py2 = 0.3*domain_radius;
elseif set == 2
     px1 = 0.1*domain_radius;
     py1 = 0;
     px2 = 0.1*domain_radius;
     py2 = 0.2*domain_radius;
 elseif set == 3
     px1 = 0;
     py1 = 0.98*domain_radius;
     px2 = 0;
     py2 = -0.98*domain_radius;
end

%%%%%%%%%%%% find random positions in the sample %%%%%%%%%%%%
if index_random==1
    for nn = 1:num_part
        while 1
        x = domain_radius*2*(0.5-rand);
        y = domain_radius*2*(0.5-rand);

        if  x^2+y^2 < domain_radius^2
            break
        end

        end
        posx_part1(nn) = x;
        posy_part1(nn) = y;
    end
else
     posx_part1(1) = px1;
     posy_part1(1) = py1;
     
     posx_part1(2) = px2;
     posy_part1(2) = py2;
end


%%%%%%%%%%%%%%%%%%%%%% PARAMETERS FORCES %%%%%%%%%%%%%%%%%%%%%%%

%parameters for the potential
%the width ot the potential/force (it is a decreasing exponential)
f_w_part        = 12;       %width potential for particle interactions
f_w_bound       = 15;       %width potential for wall interactions

%amplitude of forces for the particule and the boundaries
f0              = 0.00;    % amplitude of force for particule/particule
f_edge          = 0.007;     % amplitude of force for particule/edge

%%%%%%%%parameters for the hill function near edges
index_hill      = 1;        % bolean depending on weither you want to use the hill correction
hill_width      = 18;       % width of the hill zone
hill_power      = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETERS DYNAMICS

zeta            = 5;        % friction coeff
f_noise         = 0.0005;   % noise in the force (this is the amplitude, the direction is random at each time step/like brwonian motion)

%%%%%%%%%%%%%%%%%%%% PARAMETERS DIVISION %%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%% PARAMETERS SIMULATION %%%%%%%%%%%%%%%%%%%%%%%

dt              = 0.2;      % time step simulation
time_tot        = 300000;    % total time simulation
incr            = time_tot/dt; % total  simulation interations number
dtime_save      = 50;       % save every 1
step_save       = dtime_save/dt; % data are save every incr_save

%%%%%%%%% LOOP ON TIME %%%%%%%%
tsave           = 1;
div_index       = 0;

for tt = 1:incr
    if(mod(tt,25000)==0);disp(['on time ', num2str(tt), ' of ', num2str(incr)]); end

    time = tt*dt;

    force_calculation
    %debug if not a number
    if isnan(posx_part1(1))==1
    return
    end

end  
load('ExperCount.mat')
posx{1} = posx{2}; posy{1} = posy{2};
xx      = cat(1,posx{:});
yy      = cat(1,posy{:});

save(sprintf('./MatFiles/SingleAster_Track1_run%02d.mat',index_track),'xx','yy')        
index_track = index_track+1;
save('ExperCount.mat','index_track')

end
%% Figure generation for n  = 1
sim_iter    = 5;
x_sim       = linspace(0,1,101);
mean_finder = zeros(sim_iter,length(x_sim)-1);
figure, hold on
for i = 1:sim_iter
    load(sprintf('./MatFiles/SingleAster_Track1_run%02d.mat',i))
    rr      = sqrt(xx.^2+yy.^2);
    veloc   = abs((rr(3:end) - rr(1:end-2))/(2*dt));
    [a,b]   = data_binning1((domain_radius - rr(2:end-1))/max(domain_radius- rr(2:end-1)),veloc,x_sim);
    mean_finder(i,:)    = a(1,:);
    plot((domain_radius - rr(2:end-1))/max(domain_radius- rr(2:end-1)),veloc,'-','Color',[0.75,0.75,0.75],'LineWidth',1)
end
plot(b,mean(mean_finder/1,'omitnan'),'k-','LineWidth',3)

load('One_Aster_Profile.mat')
errorbar(Norm_x,mean(vel_x,'omitnan'),std(vel_x,'omitnan'),'ko','LineWidth',1,'MarkerSize',12)
% plot(Norm_x,mean(vel_x,'omitnan')+std(vel_x,'omitnan'),'k-','LineWidth',2)
% plot(Norm_x,mean(vel_x,'omitnan')-std(vel_x,'omitnan'),'k-','LineWidth',2)
hold off
axis([0 1 0 0.1])

%% plot time course
figure, hold on
for i = 1:sim_iter
    load(sprintf('./MatFiles/SingleAster_Track1_run%02d.mat',i))
    rr      = sqrt(xx.^2+yy.^2);
    veloc   = abs((rr(3:end) - rr(1:end-2))/(2*dt));
    timer   = dt*(1:length(veloc));
    plot(timer, veloc)
end
hold off

%% Map of movement
%%%%%%%%%%%%%%%%%% Generate figures %%%%%%%%%%%%%%%%%%
mp_time = jet(tsave);
figure, hold on
circle(0,0,domain_radius)

for tt = 1:10:tsave-1
for nn = 1:num_part

   scatter(posx{tt}(nn),posy{tt}(nn),50,mp_time(tt,:),'filled');
   hold on
   %scatter(posx_virt{tt}(nn),posy_virt{tt}(nn),50,mp_time(tt,:),'filled');
end
end

axis equal


