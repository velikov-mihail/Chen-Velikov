clear
clc

if ~exist([pwd,'Figures'], 'dir')
    mkdir(['Figures']);
    addpath([pwd,'/Figures']);
end

%% Figure 1: Anomaly mean long/short returns

clear
clc
close all

load organized_results

hold on

retall = [...
    sumall.base.insamp.g.rbar ...
    sumall.opt.insamp.n.rbar ...
    sumall.opt.postpub.n.rbar ...
    sumall.opt.postpub05.n.rbar ...
    ];
    
ret = 100*100*nanmean(retall);
SE  = 100*100*nanstd(retall)./sqrt(sum(~isnan(retall)));


defColors=get(0,'defaultAxesColorOrder');


hb = bar(ret,'stacked','BarWidth',0.5);
ha = gca;
ha.XTick = 1:4;
hb.FaceColor = defColors(1,:);


ha.XTickLabel = {
    'Gross In-Sample'
    'Net In-Sample'
    'Net Post-Pub'
    'Net Post-Pub & Post-2005'
    };

er = errorbar(ret,2*SE);

er.Color = 'k';
er.LineStyle = 'none';
er.LineWidth = 2.5;
er.MarkerSize =6;

hl=findobj(gcf,'type','legend');
delete(hl)


ylabel({'Mean Net Return', ...
        'Across 204 Anomalies (bps per month)'})
xlim([0.5 4.5]); ylim([0 80])
yyaxiscopy;

setplot
set(gcf,'position',[   934   372   895   440])


legendSize=30;
set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Position',[0 0 8.5 11])

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figure1.pdf','-transparent');


%% Figure 2: LF bias figure

clear
clc
close all

load effSpreadStruct_cv
load pdates
load dates

legendSize=45;
axisSize=40;
lineWidth=5.5;
lineWidthThin = 2;
lightGray=[0.66,0.66,0.66];
darkGray=[0.33,0.33,0.33];

    
tgood = dates <= 202012;
pdates = pdates(tgood);
 
% error in pct pt
egibbs = 100*nanmedian(effSpreadStruct.gibbs(tgood,:) - effSpreadStruct.hf_spreads.ave(tgood,:),2); 
ehl = 100*nanmedian(effSpreadStruct.hl(tgood,:) - effSpreadStruct.hf_spreads.ave(tgood,:),2); 
echl = 100*nanmedian(effSpreadStruct.chl(tgood,:) - effSpreadStruct.hf_spreads.ave(tgood,:),2); 
evov = 100*nanmedian(effSpreadStruct.vov(tgood,:) - effSpreadStruct.hf_spreads.ave(tgood,:),2); 

clf; hold on;
xlim([1983 2021])
ylim([-2 3])
line([1980 2021],[0 0],'color',darkGray);

h(1)=plot(pdates,egibbs);
h(2)=plot(pdates,ehl,'--');
h(3)=plot(pdates,echl,':');
h(4)=plot(pdates,evov,'-.');

set(h,'LineWidth',lineWidth);

defColors=get(0,'defaultAxesColorOrder');
for i=1:4
    set(h(i),'Color',defColors(i,:));
end


line([1993 1993],ylim,'color',lightGray,'lineWidth',2);
line([2003+8/12 2003+8/12], ylim,'color',lightGray,'lineWidth',2);   


hleg=legend('','Gibbs','HL','CHL','VoV','Location','Northwest');
hleg.FontSize=legendSize;
legend boxoff
ylabel({'Median (LF Spread - TAQ Spread)','(% Point)'});

yyaxiscopy


% label regions
tempy =2.3;
text(1987,tempy,'\leftarrow ISSM','FontSize',legendSize)
text(1997,tempy,'MTAQ','FontSize',legendSize)
text(2004.2,tempy,'DTAQ \rightarrow','FontSize',legendSize)

set(gca,'FontSize',axisSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figure2.pdf','-transparent');

%% Figure 3: Combined effective spreads over time

clear
clc
close all

load dates
load pdates
load tcosts_cv

tgood = dates <= 20202012;
dates = dates(tgood);
pdates = pdates(tgood);
tcosts = tcosts(tgood,:);

% find quartiles
quarts=[prctile(tcosts,[25 50 75],2)];

legendSize=45;
axisSize=30;
lineWidth=6.5;
lineWidthThin = 2;

lightGray=[0.66,0.66,0.66];

clf; hold on;
line([1983 1983],[0 12],'color',lightGray,'lineWidth',2);
line([1993 1993],[0 12],'color',lightGray,'lineWidth',2);
line([2003+8/12 2003+8/12], [0 12],'color',lightGray,'lineWidth',2);   
[H]=plot(pdates,200*quarts); % composite

set(H,'LineWidth',lineWidth);

defColors=get(0,'defaultAxesColorOrder');

% pretty up composite line
set(H(1),'LineStyle','-.','color',defColors(2,:))
set(H(2),'LineStyle','-','color',defColors(1,:))
set(H(3),'LineStyle',':','color',defColors(2,:))
AX=gca;
hleg = legend(flipud(H),{'75^{th} percentile','50^{th} percentile','25^{th} percentile'}...
    ,'location','northeast');
hleg.Position = [ 0.1573    0.7492    0.2125    0.1476];
legend boxoff
set(AX,'XLim',[pdates(1),pdates(end)]);
set(AX,'FontSize',axisSize);
ylabel(AX,'Effective spread (%)','FontSize',axisSize,'color','k')
set(gca,'YTick',0:2:20);

% label regions
tempy =7.9;
text(1978-34,tempy,'\leftarrow Average of LF proxies','FontSize',legendSize)
text(1982+1.7,tempy,'ISSM','FontSize',legendSize)
text(1993+0.4,tempy,'MTAQ','FontSize',legendSize)
text(2004+0.5,tempy,'DTAQ \rightarrow','FontSize',legendSize)

xlim([1925 2021])
ylim([0 12])
    
% label right y axis
yyaxiscopy

set(gcf,'position', [0 0  960         420])

set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figure3.pdf','-transparent');



%% Figure 4: Distribution of Spreads Paid by Academic Implementation in 2014

clear
clc
close all

load tcosts_cv
load NYSE
load universe
load dates

s=find(dates==201401);
e=find(dates==201412);
n=(e-s+1)*size(tcosts,2);

% % This is used to create tcosts_anoms
[anoms, labels]=getAnomalySignals('signed_predictors_dl_wide.csv','permno','yyyymm');
load ret
load prc
load me
anoms(:,:,203)=-abs(prc);
anoms(:,:,204)=-ret;
anoms(:,:,205)=-me;
labels(203)={'Price'};
labels(204)={'STreversal'};
labels(205)={'Size'};
stayIndicator=~strcmp(labels,'IO_ShortInterest');
labels=labels(1,stayIndicator);
anoms=anoms(:,:,stayIndicator);
tcosts_anoms=[];
for i=1:length(labels)
    fprintf('\n\n\nNow working on %s, which is %d out of %d.\n\n\n',char(labels(i)),i,length(labels));
    var=anoms(:,:,i);
    ind=makeUnivSortInd(var,5);
    ind=1*(ind==1)+2*(ind==5);
    temp=tcosts;
    temp(ind~=1 & ind~=2)=nan;
    temp=reshape(temp(s:e,:),n,1);
    temp(isnan(temp))=[];
    tcosts_anoms=[tcosts_anoms;temp];
    clear temp ind var
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(i)),i,length(labels));
end
save Results\tcosts_anoms tcosts_anoms
load tcosts_anoms

% Do all stocks
temp=tcosts;
temp=temp(s:e,:);
temp=reshape(temp,n,1);
tcosts_all=temp(isfinite(temp));

% NYSE
temp=tcosts;
temp(NYSE~=1)=nan; % Exclude non-NYSE
temp=temp(s:e,:);
temp=reshape(temp,n,1);
tcosts_nyse=temp(isfinite(temp));


% Russell 2000
temp=tcosts;
temp(universe(2).ind~=1)=nan; % Exclude non-Russell-2000
temp=temp(s:e,:);
temp=reshape(temp,n,1);
tcosts_russell=temp(isfinite(temp));

bin=100*[0:0.0001:0.045];
clf; hold on;

[f,bin] = hist(200*tcosts_anoms,bin);
hanom = plot(bin,f/sum(f));

[f,bin] = hist(200*tcosts_nyse,bin);
hnyse = plot(bin,f/sum(f));

[f,bin] = hist(200*tcosts_all,bin);
hall = plot(bin,f/sum(f));

[f,bin] = hist(200*tcosts_russell,bin);
hrus = plot(bin,f/sum(f));

xlim([0 0.8]);
lineWidth=5;
legendSize=40;

[lgd]=legend([hanom,hall,hnyse,hrus]...
    ,{'Paid by Anomaly Portfolios','All Stocks','NYSE','Russell 2000'});
legend boxoff
lgd.FontSize=legendSize;
lgd.ItemTokenSize=[60 30]';

xlabel('Effective Spread (%)'); ylabel('Frequency');


% make pretty
xlim([0 .80]);

setplot

defColors=get(0,'defaultAxesColorOrder');

set(hanom,'linewidth',2*lineWidth,'color',[ 0.5843 0.8157 0.9882])
set(hall,'linewidth',lineWidth,'linestyle','--','color', 	[0.6350, 0.0780, 0.1840])
set(hnyse,'linewidth',lineWidth,'linestyle','-.','color',defColors(3,:))
set(hrus,'linewidth',lineWidth,'linestyle',':','color',defColors(4,:))

set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figure4.pdf','-transparent')

%% Figure 5: Event-time net returns for cost-optimized implementation

clear
clc
close all

load organized_results

temp = resall.opt;

tspine = (-20*12:1:40*12)';

relret = nan(length(tspine), length(temp));

for signali = 1:length(temp)
    nret = 100*temp(signali).netpret;
    months_since_pub = temp(signali).months_since_pub;
    
    relmin = max( months_since_pub(1), tspine(1));
    relmax = min( months_since_pub(end), tspine(end));    
        
    relret(...
        tspine >= relmin & tspine <= relmax, signali ...
        ) = ...
        nret( ...
        months_since_pub >= relmin & months_since_pub <= relmax ...
        );    
end

% --- user --- %
Tsmooth = 12*5;
time        = tspine/12;

% focus on years with many anomalies
xlimnum = [-10 20]; 

xtext = 'Years Since Publication';
xticknum = xlimnum(1):5:xlimnum(2);
yticknum = -0.1:0.1:0.7;
ylimnum = [min(yticknum) max(yticknum)];


% --- plot --- %    
clf; hold on;    
      
[lez,lsmooth,lse1] = plotretbytime(time,relret,Tsmooth);    

% make pretty
xlim(xlimnum); ylim(ylimnum)        
setplot('',xtext,{'Mean Net Return','(% Monthly)'},1.25)
set(gca,'XTick',xticknum);
set(gca,'YTick',yticknum);    
yyaxis right
set(gca,'YColor','k')
ylim(ylimnum)
set(gca,'YTick',yticknum);

vline(0,'k-')

legendSize=40;
axisSize=25;
lineWidth=5;

lgd=legend([lez lsmooth lse1]...
    ,{'Average in Month','Trailing 5-year Ave','2 SE C.I.'}...
    ,'location','northeast' ...
    ,'fontsize',legendSize);
lgd.EdgeColor='w';
set(lsmooth ,'linewidth',lineWidth)
set(gcf,'position', [1000 400 560         420])

set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figure5.pdf','-transparent')

%% Figure 6: Post-publication Performance Decay Net of Costs and Stale Data

clear
clc
close all

load organized_results

ylimnum = [-100 200];
xlimnum = [-100 200];

defColors=get(0,'defaultAxesColorOrder');

xtext = -90; ytext = 100;
clf; hold on;

subplot(2,2,1); hold on;
ylim(ylimnum); xlim(xlimnum)
hline(0,'k'); vline(0,'k')
x = 1e4*sumall.base.insamp.g.rbar; y =  1e4*sumall.base.postpub.g.rbar;
model= fitlm(x,y,'intercept',false);
xfit = -1e4:1e4; yfit = model.Coefficients.Estimate*xfit;
plot(x ,y ,'x','color',defColors(1,:))
plot(xfit, yfit, '-','color',defColors(2,:))
text(xtext, ytext, {['slope = ' num2str(model.Coefficients.Estimate,2)]...
    ,['SE = ' num2str(model.Coefficients.SE,1)]}, 'color',defColors(2,:))
setplot('A. Original Gross Decay'...
    ,'In-Sample Return (bps monthly)'...
    ,'Post-Pub Return (bp p.m.)')
set(gca,'FontName','Times New Roman')


subplot(2,2,2); hold on;
ylim(ylimnum); xlim(xlimnum)
hline(0,'k'); vline(0,'k')
x = 1e4*sumall.base.insamp.g.rbar; y =  1e4*sumall.base.postpub05.g.rbar;
model= fitlm(x,y,'intercept',false);
xfit = -1e4:1e4; yfit = model.Coefficients.Estimate*xfit;
plot(x ,y ,'x','color',defColors(1,:))
plot(xfit, yfit, '-','color',defColors(2,:))
text(xtext, ytext, {['slope = ' num2str(model.Coefficients.Estimate,2)]...
    ,['SE = ' num2str(model.Coefficients.SE,1)]}, 'color',defColors(2,:))
setplot('B. Original Gross Post-2005 Decay'...
    ,'In-Sample Return (bps monthly)'...
    ,'Post-Pub05 Return (bp p.m.)')
set(gca,'FontName','Times New Roman')

subplot(2,2,3); hold on;
ylim(ylimnum); xlim(xlimnum)
hline(0,'k'); vline(0,'k')
x = 1e4*sumall.base.insamp.n.rbar; y =  1e4*sumall.base.postpub05.n.rbar;
model= fitlm(x,y,'intercept',false);
xfit = -1e4:1e4; yfit = model.Coefficients.Estimate*xfit;
plot(x ,y ,'x','color',defColors(1,:))
plot(xfit, yfit, '-','color',defColors(2,:))
text(xtext, ytext, {['slope = ' num2str(model.Coefficients.Estimate,2)]...
    ,['SE = ' num2str(model.Coefficients.SE,1)]}, 'color',defColors(2,:))
setplot('C. Original Net Post-2005 Decay'...
    ,'In-Sample Return (bps monthly)'...
    ,'Post-Pub05 Return (bp p.m.)')
set(gca,'FontName','Times New Roman')


subplot(2,2,4); hold on;
ylim(ylimnum); xlim(xlimnum)
hline(0,'k'); vline(0,'k')
x = 1e4*sumall.opt.insamp.n.rbar; y =  1e4*sumall.opt.postpub05.n.rbar;
model= fitlm(x,y,'intercept',false);
xfit = -1e4:1e4; yfit = model.Coefficients.Estimate*xfit;
plot(x ,y ,'x','color',defColors(1,:))
plot(xfit, yfit, '-','color',defColors(2,:))
text(xtext, ytext, {['slope = ' num2str(model.Coefficients.Estimate,1)]...
    ,['SE = ' num2str(model.Coefficients.SE,1)]}, 'color',defColors(2,:))
setplot('D. Cost-Mitigated Net Post-2005 Decay'...
    ,'In-Sample Return (bps monthly)'...
    ,'Post-Pub05 Return (bp p.m.)')

set(gcf,'position', [1000 400 775         608])

set(gca,'FontName','Times New Roman')
export_fig('Figures/figure6.pdf','-transparent')

%% Figure 7: Heterogeneity in cost-optimized mean net return post-publication and post-2005

clear
clc
close all

load organized_results
retsamp = sumall.opt.postpub05.n.rbar*100*100;
tstatsamp  = sumall.opt.postpub05.n.tstat;
namesamp = anomaly_summary.Acronym;

% clean up names
for i = 1:length(namesamp)
    % remove underscores    
    namesamp{i}(namesamp{i} == '_') = [];    
    % crop
    namesamp{i} = namesamp{i}(1:min(length(namesamp{i}),8));
end



% choose qcut = 0.25;
tempred = 'red';
temppurple = [1 0 1]*0.6; % purple
tempblue = 'blue';

clf; hold on;

strmax = 9;
fontsize = 12.5;

bin = -100:20:120;

binloc = 1:(length(bin)-1);

% produce frame
ymax = ceil((max(histcounts(retsamp,bin))+1)/10)*10;

xlim([0.8 binloc(end)+1.1]);  ylim([0 ymax])

% plot
for bini = 1:length(bin)-1
    
    jlist = find(retsamp >= bin(bini) & retsamp < bin(bini+1));
    
    % make alphabetical
    [namealpha,ialpha] = sort(namesamp(jlist));    
    namealpha = flipud(namealpha);
    ialpha = flipud(ialpha);
    
    for portj = 1:length(jlist)  
        % crop string
        str = namealpha{portj};
        tempn = min(length(str),strmax); 
        str = str(1:tempn);
        
        % highlight subsets (choose color)
        tempcolor = 'black';
        tempangle = 'normal';
        tempweight = 'normal';
        
        % convert to standard index
        porti = jlist(ialpha(portj));        

        if abs(tstatsamp(porti)) > 2
            tempcolor = tempblue;
            tempweight = 'bold';     
            str = ['\textbf{' str '}'];
        end      
        
        if abs(tstatsamp(porti)) < 1.5
            tempcolor = tempred;
            tempangle = 'italic';     
            str = ['\emph{' str '}'];
        end        
        
        % write down name in hist
        text(binloc(bini), ...
            portj,...
            str,...
            'color',tempcolor,...
            'fontsize',fontsize,...
            'fontangle',tempangle,...
            'fontweight',tempweight,...
            'interpreter','latex');    
    end    
end

% create xticks
h = get(gca);
h.XAxis.TickLength = [0 0]; % remove ticks
h.XAxis.TickLabel = bin(1:end);

% label
setplot('','Net Return Post-Pub and Post-2005 (bps per month)','Count',1.3,'helvetica','none')

% vertical line
bin0 = find(bin == 0) - 0.03;
plot([bin0 bin0],[0 60],'k:','linewidth', 1.5);
h.XAxis.TickValues = [binloc binloc(end)+1]+0.05; % reposition xticks 


set(gca,'FontSize',25);

% legend alt
tempx = 1.5;
tempy = 45;
legsize = 20;
text(tempx,tempy-0*1.4,'\textit{italics: t-stat $<$ 1.5}' ...
    ,'interpreter','latex','color',tempred,'fontsize',legsize)
text(tempx,tempy-2*1.4,'\textbf{bold: t-stat $>$ 2.0}' ...
    ,'interpreter','latex','color',tempblue,'fontsize',legsize)

set(gcf,'position', [      965   200   1500   1100])

set(gca,'FontName','Times New Roman')
export_fig('Figures/figure7.pdf','-transparent')


%% Figure 8: Dollar invested in combo strategies

clear
clc
close all

load res_combo_strats

nstrats=size(resISew,2);
ncombos=size(resISew,1);

netxretIS=nan(ncombos,nstrats);
netxretOS=nan(ncombos,nstrats);


for i=1:size(resOSew,1)
    for j=1:size(resOSew,2)
        netxretIS(i,j)=resISew(i,j).netxret; 
        netxretIS(i,j+nstrats)=resISvw(i,j).netxret;
        netxretOS(i,j)=resOSew(i,j).netxret; 
        netxretOS(i,j+nstrats)=resOSvw(i,j).netxret;                       
    end
end

for r=1:4
    cr(r,1)=find(netxretIS(r,:)==max(netxretIS(r,:)));
end

resew=struct;

for i=1:4
    res(i,1).dates=[resISew(i,1).dates(:,end); resOSew(i,1).dates(:,end)];
    res(i,1).pretew=[resISew(i,1).pret(:,end); resOSew(i,1).pret(:,end)];
    res(i,1).pretvw=[resISvw(i,1).pret(:,end); resOSvw(i,1).pret(:,end)];    
    res(i,1).pretew_opt=[resISew(i,cr(i)).pret(:,end); resOSew(i,cr(i)).pret(:,end)];
    res(i,1).pretvw_opt=[resISvw(i,cr(i)).pret(:,end); resOSvw(i,cr(i)).pret(:,end)];    
    res(i,1).netpretew=[resISew(i,1).netpret(:,end); resOSew(i,1).netpret(:,end)];
    res(i,1).netpretvw=[resISvw(i,1).netpret(:,end); resOSvw(i,1).netpret(:,end)];    
    res(i,1).netpretew_opt=[resISew(i,cr(i)).netpret(:,end); resOSew(i,cr(i)).netpret(:,end)];
    res(i,1).netpretvw_opt=[resISvw(i,cr(i)).netpret(:,end); resOSvw(i,cr(i)).netpret(:,end)];    
end

dates=res(1).dates;
pdates = floor(dates/100) + mod(dates,100)/12;


fontSize=15;
lineWidth=3.5;

% Check if daily or monthly returns
if length(num2str(dates(1))) == 8
    x = floor(dates/10000) + (floor(mod(dates,10000)/100)-1)/12 + mod(dates,100)/(12*30);
else
    x = floor(dates/100) + mod(dates,100)/12;
end

load OOSreturns_Fig6
dates2=floor(dates2/100);
dnmu=nan(size(dates));
[~,ia,ib]=intersect(dates,dates2);
dnmu(ia)=ReturnsOOS2(ib,1)-ReturnsOOS2(ib,4);


figure;

subplot(2,1,2);
rets=[res.netpretew_opt dnmu];
x = [2*x(1)-x(2); x];
rets=rets.*repmat(nanstd(rets,[],1),size(rets,1),1)/nanstd(rets(:,1));
rets(isnan(rets)) = 0;
y = cumprod(1+[zeros(1,size(rets,2)); rets]);
b=(find(dates==dates2(1))-1);

defColors=get(0,'defaultAxesColorOrder');

h=semilogy(x,y(:,1),'-',x,y(:,2),'--',x,y(:,3),':',x,y(:,4),'-','LineWidth',lineWidth);
for i=1:4
    set(h(i),'Color',defColors(i,:));
end
xlim([x(1) x(end)])
ylim([0 200]);
set(gca,'YTick',10.^[0:8])
set(gca,'YTickLabel',{'$1';'$10';'$100';'$1,000';'$10,000';'$100,000';'$1,000,000';'$100,000,000'})
p=xline(pdates(dates==200512),'lineWidth',2);
uistack(p,'bottom');
title({'B. With cost-mitigation optimized over 1985-2005'},'FontWeight','normal');
lgd=legend(h(1:4),'Fama-MacBeth','IPCA','Average rank','LASSO','Location','Northwest');
legend boxoff
lgd.FontSize=fontSize;
set(gca,'FontSize',fontSize);
set(gca,'FontName','Times New Roman')


subplot(2,1,1);
rets=[res.netpretew dnmu];
rets=rets.*repmat(nanstd(rets,[],1),size(rets,1),1)/nanstd(rets(:,1));
rets(isnan(rets)) = 0;
y = cumprod(1+[zeros(1,size(rets,2)); rets]);
b=(find(dates==dates2(1))-1);
y(1:b,5)=nan;
y((find(dates==dates2(end))+1):end,5)=nan;
y(:,5)=y(:,5)+mean(y(b+1,1:4))-1;
h=semilogy(x,y(:,1),'-',x,y(:,2),'--',x,y(:,3),':',x,y(:,4),'-.',x,y(:,5),'-','LineWidth',lineWidth);
set(h(5),'Color',[0.8 0.8 0.8]);
uistack(h(5),'bottom');
xlim([x(1) x(end)])
ylim([0 200]);
set(gca,'YTick',10.^[0:8])
set(gca,'YTickLabel',{'$1';'$10';'$100';'$1,000';'$10,000';'$100,000';'$1,000,000';'$100,000,000'})
p=xline(pdates(dates==200512),'lineWidth',2);
uistack(p,'bottom');
title({'A. Net equal-weighted decile sorts on fitted expected returns'},'FontWeight','normal');
set(gca,'FontSize',fontSize);
lgd=legend(h(1:5),'Fama-MacBeth','IPCA','Average rank','LASSO','DMNU','Location','Northwest');
legend boxoff
lgd.FontSize=fontSize;
set(gca,'FontName','Times New Roman')
orient(gcf,'landscape')

set(gcf,'position', [      965   200   800   1207])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

export_fig('Figures/figure8.pdf','-transparent');

%% Figure A.1: Distribution of Net Returns: In-Sample, Before Cost Optimization

clear
clc
close all

load organized_results

% choose varaibles
temp = sumall.base;
retsamp = temp.insamp.n.rbar*100*100;
boldme = temp.insamp.basic.turnover*100 > 30;
emphme = zeros(size(retsamp));
namesamp = anomaly_summary.Acronym;

% clean up names
for i = 1:length(namesamp)
    % remove underscores    
    namesamp{i}(namesamp{i} == '_') = [];    
    % crop
    namesamp{i} = namesamp{i}(1:min(length(namesamp{i}),8));
end


qcut = 0.25;
tempred = 'red';
temppurple = [1 0 1]*0.6; % purple
tempblue = 'blue';

clf; hold on;

strmax = 9;
fontsize = 12.5;

bin = [-500 -300 -100 -80:20:100 300];

binloc = 1:(length(bin)-1);

% produce frame
ymax = ceil((max(histcounts(retsamp,bin))+1)/10)*10;

xlim([0.8 binloc(end)+1.1]);  ylim([0 ymax])

% plot
for bini = 1:length(bin)-1
    
    jlist = find(retsamp >= bin(bini) & retsamp < bin(bini+1));
    
    % make alphabetical
    [namealpha,ialpha] = sort(namesamp(jlist));    
    namealpha = flipud(namealpha);
    ialpha = flipud(ialpha);
    
    for portj = 1:length(jlist)  
        % crop string
        str = namealpha{portj};
        tempn = min(length(str),strmax); 
        str = str(1:tempn);
        
        % highlight subsets (choose color)
        tempcolor = 'black';
        tempangle = 'normal';
        tempweight = 'normal';
        
        % convert to standard index
        porti = jlist(ialpha(portj));        

        if boldme(porti)
            tempcolor = tempblue;
            tempweight = 'bold';     
            str = ['\textbf{' str '}'];
        end      
        
        if emphme(porti)
            tempcolor = tempred;
            tempangle = 'italic';     
            str = ['\emph{' str '}'];
        end        
        
        
        % write down name in hist
        text(binloc(bini), ...
            portj,...
            str,...
            'color',tempcolor,...
            'fontsize',fontsize,...
            'fontangle',tempangle,...
            'fontweight',tempweight,...
            'interpreter','latex');    
    end    
end


% create xticks
h = get(gca);
h.XAxis.TickLength = [0 0]; % remove ticks
h.XAxis.TickLabel = bin(1:end);

% label
setplot('','Net Return In-Sample (bps per month)','Count',1.3,'helvetica','none')

% vertical line
h.XAxis.TickValues = [binloc binloc(end)+1]-0.15; % reposition xticks 


text(find(bin==-300)-0.05, 0, '//')
text(find(bin==-100)-0.05, 0, '//')
text(find(bin==100)-0.05, 0, '//')

set(gca,'FontSize',25);

% legend alt
tempx = 1.5;
tempy = 29;
legsize = 20;

text(tempx,tempy-1*1.4,'\textbf{bold: turnover $>$ 30\%}' ...
    ,'interpreter','latex','color',tempblue,'fontsize',legsize)

set(gcf,'position', [      965   200   988   707])

set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figureA1.pdf','-transparent')

%% Figure A.2: Distribution of Net Returns: In-Sample, After Cost Optimization

clear
clc
close all

load organized_results

temp = sumall.opt;

retsamp = temp.insamp.n.rbar*100*100;
boldme = strcmp(temp.set.sweight,'vw');
namesamp = anomaly_summary.Acronym;

% clean up names
for i = 1:length(namesamp)
    % remove underscores    
    namesamp{i}(namesamp{i} == '_') = [];    
    % crop
    namesamp{i} = namesamp{i}(1:min(length(namesamp{i}),8));
end

qcut = 0.25;
tempred = 'red';
temppurple = [1 0 1]*0.6; % purple
tempblue = 'blue';

clf; hold on;

strmax = 9;
fontsize = 12.5;
legfont = 20;

bin = [-40:20:160 300];

binloc = 1:(length(bin)-1);

% produce frame
ymax = ceil((max(histcounts(retsamp,bin))+1)/10)*10;

xlim([0.8 binloc(end)+1.1]);  ylim([0 ymax])

% plot
for bini = 1:length(bin)-1
    
    jlist = find(retsamp >= bin(bini) & retsamp < bin(bini+1));
    
    % make alphabetical
    [namealpha,ialpha] = sort(namesamp(jlist));    
    namealpha = flipud(namealpha);
    ialpha = flipud(ialpha);
    
    for portj = 1:length(jlist)  
        % crop string
        str = namealpha{portj};
        tempn = min(length(str),strmax); 
        str = str(1:tempn);
        
        % highlight subsets (choose color)
        tempcolor = 'black';
        tempangle = 'normal';
        tempweight = 'normal';
        
        % convert to standard index
        porti = jlist(ialpha(portj));        

        if boldme(porti)
            tempcolor = tempblue;
            tempweight = 'bold';     
            str = ['\textbf{' str '}'];
        end              
        
        % write down name in hist
        text(binloc(bini), ...
            portj,...
            str,...
            'color',tempcolor,...
            'fontsize',fontsize,...
            'fontangle',tempangle,...
            'fontweight',tempweight,...
            'interpreter','latex');    
    end    
end

% create xticks
h = get(gca);
h.XAxis.TickLength = [0 0]; % remove ticks
h.XAxis.TickLabel = bin(1:end);

% label
setplot('','Net Return In-Sample (bps per month)','Count',1.3,'helvetica','none')

h.XAxis.TickValues = [binloc binloc(end)+1]-0.15; % reposition xticks 

set(gca,'FontSize',25);

% legend alt
tempx = 6.5;
tempy = 60;
legsize = 20;

text(tempx,tempy-1*1.4,'\textbf{bold: value-weighted}' ...
    ,'interpreter','latex','color',tempblue,'fontsize',legsize)

set(gcf,'position', [      965   200   988   707])

set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))

orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  
export_fig('Figures/figureA2.pdf','-transparent')

%% Figure A.3: Distribution of Publication Years

clear
clc
close all

load organized_results

year = floor(sumall.base.set.pubdate/100);
T = sumall.base.postpub.basic.T;

clf; 
subplot(2,1,1)
histogram(year)
xlabel('Publication Year')
ylabel('# of Anomalies')
set(gca,'FontSize',25);

subplot(2,1,2)
histogram(T,30)
xlabel('Post-Pub Sample Length (Months)')
ylabel('# of Anomalies')


set(gcf,'position', [ 838   199   722   665])

set(gca,'FontSize',25);
set(gca,'FontName','Times New Roman')

 print(gcf, '-dpdf', 'Figures/figureA3.pdf','-fillpage'); 


close all