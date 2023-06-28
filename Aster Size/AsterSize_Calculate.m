% For analysing aster size

load('AsterDetails.mat')
figure,hold on
plot(b,nanmean(astersize),'bo','MarkerSize',10)
plot(b,nanmean(astersize)+nanstd(astersize,[],1),'b--','Linewidth',2)
plot(b,nanmean(astersize)-nanstd(astersize,[],1),'b--','Linewidth',2)
fitter = fit(b',nanmean(astersize)','a*exp(-(x-0.5)/b)+c','StartPoint',[1,10,0.05])
% plot(b,0.9*exp(-b/8)+0.05,'k-')
plot(b,(fitter.a)*exp(-(b-0.5)/fitter.b)+fitter.c,'k-','Linewidth',2)
hold off
ylabel('Normalised MT signal')
xlabel('Distance (\mu{m})')
set(gca,'FontSize',14)
%%

figure,plot(diff(nanmean(astersize)))

%%
range   = [3, length(b)];
asty    = nanmean(astersize);
dasty   = nanstd(astersize,[],1);
asty    = asty(range(1):1:range(2));
dasty   = dasty(range(1):1:range(2));
xx      = b(range(1):1:range(2));
figure,hold on
plot(b,nanmean(astersize),'ro','MarkerSize',10)
plot(xx,asty,'bo','MarkerSize',10)
plot(xx,asty+dasty,'b--','Linewidth',2)
plot(xx,asty-dasty,'b--','Linewidth',2)
fitter = fit(xx',asty','a*exp(-(x-0.5)/b)','StartPoint',[1,10])
% plot(b,0.9*exp(-b/8)+0.05,'k-')
plot(b,(fitter.a)*exp(-(b-0.5)/fitter.b),'k-','Linewidth',2)
hold off
ylabel('Normalised MT signal')
xlabel('Distance (\mu{m})')
set(gca,'FontSize',14)

