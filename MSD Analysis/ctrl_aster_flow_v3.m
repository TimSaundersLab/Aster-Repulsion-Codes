%% Analyse no aster vs aster data

% Data from single aster-edge repulsion
dx_all      = [0.1083333,0.1083333,0.1083333,0.1083333,0.1083333,0.1625,0.1083333];
dt_all      = [30,30,15,15,15,15,15];
FigPlot     = 0;
TrackStore  = struct([]);
Window      = [1,2,3,4,6,8,12,16,24,32,48,64];
Window1     = [1:4,6,8,12,16,24,32];
MSDLength   = 128;

for Droplet = 1:7
    dx      = dx_all(Droplet);                      % Space step
    dt      = dt_all(Droplet);                      % Time step
    load(sprintf('./MatFiles/Drop%01dTracks.mat',Droplet));    % Load tracks
    Tracks      = Tracker(end,1)+1;
    for i = 1:Tracks
        tmp     = Tracker(:,1) == (i-1);
        TrackStore{Droplet}(i).Track = Tracker(tmp,3:4);
        TrackStore{Droplet}(i).Time  = Tracker(tmp,2)*dt_all*(Droplet);
        tmp     = TrackStore{Droplet}(i).Track;
        msd     = zeros(length(Window),1);
        if(length(tmp)>20)
        if(Droplet<=2)
            for k = 1:length(Window1)
                if(length(tmp)>=Window1(k))
                    if(k==1)
                        msd(2) = mean((TrackStore{Droplet}(i).Track(1:end-Window1(k),1)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,1)).^2 + ...
                           (TrackStore{Droplet}(i).Track(1:end-Window1(k),2)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,2)).^2);
                    elseif(k==2)
                        msd(4) = mean((TrackStore{Droplet}(i).Track(1:end-Window1(k),1)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,1)).^2 + ...
                           (TrackStore{Droplet}(i).Track(1:end-Window1(k),2)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,2)).^2);
                    else
                        msd(k+2) = mean((TrackStore{Droplet}(i).Track(1:end-Window1(k),1)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,1)).^2 + ...
                           (TrackStore{Droplet}(i).Track(1:end-Window1(k),2)-TrackStore{Droplet}(i).Track(Window1(k)+1:end,2)).^2);
                    end
                end
            end

        else
            msd     = zeros(length(Window),1);
        
            for k = 1:length(Window)
                if(length(tmp)>=Window(k))
                    msd(k) = mean((TrackStore{Droplet}(i).Track(1:end-Window(k),1)-TrackStore{Droplet}(i).Track(Window(k)+1:end,1)).^2 + ...
                           (TrackStore{Droplet}(i).Track(1:end-Window(k),2)-TrackStore{Droplet}(i).Track(Window(k)+1:end,2)).^2);
                end
            end

        end
        end
        TrackStore{Droplet}(i).msd  = msd;
    end
    TrackStore{Droplet}(1).Dx       = dx;
    TrackStore{Droplet}(1).Dt       = dt;
end

%% Data from no aster dynamics
dx_all      = [0.1083333,0.1083333,0.1083333,0.1083333,0.1083333,0.1083333];
dt_all      = [15,15,15,15,15,15];
TrackStore1 = struct([]);

for Droplet = 1:6
    load(sprintf('./MatFiles/Drop%01dTracks_NoAster.mat',Droplet));    % Load tracks
    Tracks      = Tracker(end,1)+1;
    for i = 1:Tracks
        tmp     = Tracker(:,1) == (i-1);
        TrackStore1{Droplet}(i).Track = Tracker(tmp,3:4);
        TrackStore1{Droplet}(i).Time  = Tracker(tmp,2);
        tmp     = TrackStore1{Droplet}(i).Track;
        msd     = zeros(length(Window),1);
        if(length(tmp)>20)
        for k = 1:length(Window)
            if(length(tmp)>=Window(k))
                msd(k) = mean((TrackStore1{Droplet}(i).Track(1:end-Window(k),1)-TrackStore1{Droplet}(i).Track(Window(k)+1:end,1)).^2 + ...
                       (TrackStore1{Droplet}(i).Track(1:end-Window(k),2)-TrackStore1{Droplet}(i).Track(Window(k)+1:end,2)).^2);
            end
        end
        end
        TrackStore1{Droplet}(i).msd  = msd;
    end
    TrackStore1{Droplet}(1).Dx       = dx;
    TrackStore1{Droplet}(1).Dt       = dt;
end

%% Analyse MSD distribution


msd1    = zeros(5,length(Window));
msd2    = zeros(6,length(Window));
figure, hold on
for i =3:7
    tmp1    = cat(2,TrackStore{i}.msd);
    for k = 1:length(Window)
        tmp2 = tmp1(k,tmp1(k,:)~=0);
        msd1(i-2,k) = mean(tmp2,'omitnan');       
    end
    plot(Window,msd1(i-2,:),'-','Color',[0.7,0.7,0.7])
end
for i = 2:6
    tmp3    = cat(2,TrackStore1{i}.msd);
    for k = 1:length(Window)
        tmp4 = tmp3(k,tmp3(k,:)~=0);
        msd2(i,k) = mean(tmp4,'omitnan');
    end
    plot(Window,msd2(i,:),'--','Color',[0.4,0.4,0.4])
end
[fitting,G0]     = fit(Window',mean(msd1,'omitnan')','a*x^r','StartPoint',[1,1]);
[fitting1,G1]    = fit(Window',mean(msd2,'omitnan')','a*x^r','StartPoint',[1,1]);
plot(Window,(mean(msd1,'omitnan')),'ko','MarkerSize',12,'LineWidth',2)
plot(Window,(mean(msd2,'omitnan')),'kd','MarkerSize',12,'LineWidth',2)
plot(Window,fitting.a*(Window).^fitting.r,'k-','LineWidth',3)
plot(Window,fitting1.a*(Window).^fitting1.r,'k--','LineWidth',3)
hold off
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
set(gca,'LineWidth',1)
% set(gca,'XScale','log')
% set(gca,'YScale','log')
