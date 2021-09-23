%% Hasbrouck(2009) Figure 3, Panel A Replication

clear
clc

set(0,'DefaultAxesFontName', 'Times New Roman')

titleSize=30;
axisSize=25;
legendSize=20;
lineWidth=5;
markerSize=5;

veryLightGray=[0.77,0.77,0.77];
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.3];
veryDarkGray=[0.22,0.22,0.22];

load pdates
load dates
load exchcd
load me
load gibbs_cv
load NYSE

indExchToExcl=exchcd~=2 | isnan(gibbs);
gibbs(indExchToExcl)=nan;
me(dates-100*floor(dates/100)~=12,:)=nan;
gibbs(dates-100*floor(dates/100)~=12,:)=nan;

ind=aprts(me,4,NYSE);

qrtile1=gibbs;
qrtile2=gibbs;
qrtile3=gibbs;
qrtile4=gibbs;

qrtile1(ind~=1)=nan;
qrtile2(ind~=2)=nan;
qrtile3(ind~=3)=nan;
qrtile4(ind~=4)=nan;

pdates=pdates-0.5;

y=[(nanmean(qrtile1,2)) ...
        (nanmean(qrtile2,2)) ...
        (nanmean(qrtile3,2)) ...
        (nanmean(qrtile4,2))];
x=floor(dates/100);
figure;

ind=isfinite(x) & isfinite(sum(y,2));
x=x(ind);
y=y(ind,:);
hold on;
plot(x,y(:,1),'-','color',veryLightGray,'LineWidth',lineWidth)
plot(x,y(:,2),'LineStyle','-.','color',lightGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
plot(x,y(:,3),'LineStyle','-','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
plot(x,y(:,4),'LineStyle','--','color',veryDarkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
ylim([0 nanmax(nanmax(y))]);
hold off;
set(gca,'FontSize',axisSize);
title('Hasbrouck Figure 3, Panel A Replication (My Data)','FontSize',titleSize);
ylabel('Effective spread','FontSize',axisSize);
legend('Quartile 1','Quartile 2','Quartile 3','Quartile 4');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

%% Corwin and Schultz Figure 2


clear
clc

set(0,'DefaultAxesFontName', 'Times New Roman')

titleSize=30;
axisSize=25;
legendSize=20;
lineWidth=5;
markerSize=5;

veryLightGray=[0.77,0.77,0.77];
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.3];
veryDarkGray=[0.22,0.22,0.22];

load pdates
load dates
load exchcd
load me
load NYSE 
load hl_cv

indExchToExcl=exchcd~=1 & exchcd~=2;

me(indExchToExcl)=nan;
hl(indExchToExcl)=nan;

ind=makeUnivSortInd(me,10);

small=hl;
large=hl;
small(ind~=1)=nan;
large(ind~=10)=nan;


pdates=pdates-0.5;

y=[(nanmean(hl,2)) ...
        (nanmean(small,2)) ...
        (nanmean(large,2))];
x=pdates;
% plot(x,y)
figure;
b=find(dates==200001);
e=find(dates==202012);
hold on;
plot(x,y(:,1),'-','color',veryDarkGray,'LineWidth',lineWidth)
plot(x,y(:,2),'LineStyle','-.','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
plot(x,y(:,3),'LineStyle','-','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
xlim([pdates(b) pdates(e)])
ylim([0 nanmax(nanmax(y(b:end,:)))+0.01]);
hold off;
set(gca,'FontSize',axisSize);
title('Corwin-Schultz Figure 2, Panel A Replication (My Data)','FontSize',titleSize);
ylabel('Effective spread','FontSize',axisSize);
legend('NYSE','Small Stocks','Large Stocks');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

%% Plot Figure 8 from Abdi and Ranaldo

clear
clc

set(0,'DefaultAxesFontName', 'Times New Roman')

titleSize=30;
axisSize=25;
legendSize=20;
lineWidth=5;
markerSize=5;

veryLightGray=[0.77,0.77,0.77];
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.3];
veryDarkGray=[0.22,0.22,0.22];

load pdates
load dates
load exchcd
load me
load chl_cv


indExchToExcl=exchcd~=2 | isnan(chl);
me(indExchToExcl)=nan;
chl(indExchToExcl)=nan;

ind=aprts(me,10);
small=chl;
large=chl;
small(ind~=1)=nan;
large(ind~=10)=nan;

pdates=pdates-0.5;

y=[(nanmean(chl,2)) ...
        (nanmean(small,2)) ...
        (nanmean(large,2))];
x=pdates;

b=find(dates==196301);
e=find(dates==201512);
figure;
hold on;
plot(x,y(:,1),'-','color',veryDarkGray,'LineWidth',lineWidth)
plot(x,y(:,2),'LineStyle','-.','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
plot(x,y(:,3),'LineStyle','-','color',darkGray,'LineWidth',lineWidth,'MarkerSize',markerSize)
xlim([pdates(b) pdates(e)])
% ylim([0 1000]);
set(gca,'FontSize',axisSize);
title('Abdi-Ranaldi Figure 8 Replication (My Data)','FontSize',titleSize);
ylabel('Effective spread','FontSize',axisSize);
legend('NYSE','Small Stocks','Large Stocks');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
























