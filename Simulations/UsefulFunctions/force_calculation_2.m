dist_part   = zeros(num_part,num_part);
amp_force   = zeros(num_part,num_part);
vect_force_x= zeros(num_part,num_part);
vect_force_y= zeros(num_part,num_part);
posx_circle = zeros(num_part,1);
posy_circle = zeros(num_part,1);
posx_virtual = zeros(num_part,1);
posy_virtual = zeros(num_part,1);
dist_virtual = zeros(num_part,num_part);
vect_force_virt_x = zeros(num_part,1);
vect_force_virt_y = zeros(num_part,1);
posx_part2 = zeros(num_part,1);
posy_part2 = zeros(num_part,1);
amp_force_virt = zeros(num_part,1);

for nn = 1:num_part
    for nnn = 1:num_part
        %calculate distance between particules
        dist_part(nn,nnn) = sqrt((posx_part1(nn)-posx_part1(nnn))^2+(posy_part1(nn)-posy_part1(nnn))^2);

        if nn == nnn
            %the particle can not exert a force on itself
            amp_force(nn,nnn) = 0;
            vect_force_x(nn,nnn) = 0;
            vect_force_y(nn,nnn) = 0;
        else
            %amplitude of the force between the two particles
            amp_force(nn,nnn)    = ((dist_part(nn,nnn))^hill_power/(hill_width_ast^hill_power+(dist_part(nn,nnn))^hill_power)) * f0*exp(-dist_part(nn,nnn)/f_w_part); %
            
            vect_force_x(nn,nnn) = (posx_part1(nn)-posx_part1(nnn))/dist_part(nn,nnn);
            vect_force_y(nn,nnn) = (posy_part1(nn)-posy_part1(nnn))/dist_part(nn,nnn);
        end  
    end
end


%distance matrix for virtual particles
for nn = 1:num_part

    dist = sqrt(posx_part1(nn)^2+posy_part1(nn)^2);

    posx_circle(nn) = domain_radius*posx_part1(nn)/dist;
    posy_circle(nn) = domain_radius*posy_part1(nn)/dist;

    dist_circle = sqrt((posx_circle(nn)-posx_part1(nn))^2+(posy_circle(nn)-posy_part1(nn))^2);

    posx_virtual(nn) = (domain_radius+dist_circle)*posx_part1(nn)/dist;
    posy_virtual(nn) = (domain_radius+dist_circle)*posy_part1(nn)/dist;

    dist_virtual(nn,nn) = sqrt((posx_virtual(nn)-posx_part1(nn))^2+(posy_virtual(nn)-posy_part1(nn))^2);
    vect_force_virt_x(nn) = (posx_part1(nn)-posx_virtual(nn))/dist_virtual(nn,nn);
    vect_force_virt_y(nn) = (posy_part1(nn)-posy_virtual(nn))/dist_virtual(nn,nn);

    if index_hill == 1
        amp_force_virt(nn) = (dist_virtual(nn,nn)^hill_power/(hill_width^hill_power+dist_virtual(nn,nn)^hill_power))*f_edge*exp(-dist_virtual(nn,nn)/f_w_bound);
    else
        amp_force_virt(nn) = f_edge*exp(-dist_virtual(nn,nn)/f_w_bound);
    end


end

%calculte all the forces that apply to one particle and then sum of
%forces
for nn = 1:num_part
    posx_part2(nn) = posx_part1(nn)+ dt*(sum(amp_force(nn,:).*vect_force_x(nn,:)))/zeta + dt*(amp_force_virt(nn)*vect_force_virt_x(nn))/zeta +dt*(f_noise*2*(0.5-rand))/zeta;
    posy_part2(nn) = posy_part1(nn)+ dt*(sum(amp_force(nn,:).*vect_force_y(nn,:)))/zeta + dt*(amp_force_virt(nn)*vect_force_virt_y(nn))/zeta +dt*(f_noise*2*(0.5-rand))/zeta;      
end

%     for nn = 1:num_part
%         posx_part2(nn) = posx_part1(nn)+ dt*(sum(amp_force(nn,:).*vect_force_x(nn,:)))/zeta ;
%         posy_part2(nn) = posy_part1(nn)+ dt*(sum(amp_force(nn,:).*vect_force_y(nn,:)))/zeta ;      
%     end

if tt == 1
    posx{tsave} = posx_part1;posy{tsave} = posy_part1;
    posx_virt{tsave} = posx_virtual;posy_virt{tsave} = posy_virtual;
    time_save(tsave) = time;
    tsave = tsave+1;
end
if mod(tt,step_save)==0
    posx{tsave} = posx_part1;posy{tsave} = posy_part1;
     posx_virt{tsave} = posx_virtual;posy_virt{tsave} = posy_virtual;
     time_save(tsave) = time;
    tsave = tsave+1;
end


posx_part1 = posx_part2;posy_part1 = posy_part2;