load('Four_Aster_Angles.mat')
Ang_Ran = linspace(0,180,37);

% N_Rand  = 1000;
% AngleRan1   = pi * rand(N_Rand,1);
% AngleRan2   = (pi - AngleRan1).*rand(N_Rand,1);
% AngleRan3   = pi - AngleRan1 - AngleRan2;

Angles1 = max(tmp')';

[a,~]   = histc(Angles1,Ang_Ran);
% [a1,~]  = histc(Angles_Noisy,Ang_Ran);
% [a1,~]  = histc(180*cat(1,AngleRan1,AngleRan2,AngleRan3)/pi,Ang_Ran);
figure,hold on
plot(Ang_Ran(1:end)+2.5,a/sum(a),'bo-')
% plot(Ang_Ran(1:end)+2.5,a1/sum(a1),'ro-')   % This is same parameters but with noisy 100x larger
hold off
