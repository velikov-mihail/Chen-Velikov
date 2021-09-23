%% Start a diary (will have all the table results)
clear
clc

diary Results/ChenVelikovTablesOutput.txt

fprintf('\n\n\n\nTable printing started @ %s\n\n\n',char(datetime('now')));

%% Table 1: Correlations between effective Bid-ask spreads

fprintf('\n\n\nTable 1 output:\n\n\n');

clear
clc

load effSpreadStruct_cv
load dates

n=size(effSpreadStruct.gibbs,1)*size(effSpreadStruct.gibbs,2);
reshapedSpreads=[reshape(effSpreadStruct.gibbs,n,1) ... 
    reshape(real(effSpreadStruct.hl),n,1) ... 
    reshape(effSpreadStruct.chl,n,1) ... 
    reshape(effSpreadStruct.vov,n,1)];

reshapedSpreads(isnan(sum(reshapedSpreads,2)),:)=[];
corrs=triu(corrcoef(reshapedSpreads))';
corrs(corrs==0)=nan;
mat2Tex(corrs,corrs,{'Gibbs','HL','CHL','VoV'},2);
clear reshapedSpreads corrs

s=find(dates==199301);
n=size(effSpreadStruct.gibbs(s:end,:),1)*size(effSpreadStruct.gibbs,2);
reshapedSpreads=[reshape(effSpreadStruct.hf_spreads.ave(s:end,:),n,1) ...
    reshape(effSpreadStruct.gibbs(s:end,:),n,1) ... 
    reshape(real(effSpreadStruct.hl(s:end,:)),n,1) ... 
    reshape(effSpreadStruct.chl(s:end,:),n,1) ... 
    reshape(effSpreadStruct.vov(s:end,:),n,1) ...
    ];
reshapedSpreads=[reshapedSpreads nanmean(reshapedSpreads(:,2:5),2)];


reshapedSpreads(isnan(sum(reshapedSpreads,2)),:)=[];
corrs=triu(corrcoef(reshapedSpreads))';
corrs(corrs==0)=nan;
mat2Tex(corrs,corrs,{'TAQ','Gibbs','HL','CHL','VoV','LF Ave'},2);


s=find(dates==198301);
e=find(dates==199212);
n=size(effSpreadStruct.gibbs(s:e,:),1)*size(effSpreadStruct.gibbs,2);
reshapedSpreads=[reshape(effSpreadStruct.hf_spreads.ave(s:e,:),n,1) ...
    reshape(effSpreadStruct.gibbs(s:e,:),n,1) ... 
    reshape(real(effSpreadStruct.hl(s:e,:)),n,1) ... 
    reshape(effSpreadStruct.chl(s:e,:),n,1) ... 
    reshape(effSpreadStruct.vov(s:e,:),n,1) ...
    ];

reshapedSpreads=[reshapedSpreads nanmean(reshapedSpreads(:,2:5),2)];


reshapedSpreads(isnan(sum(reshapedSpreads,2)),:)=[];
corrs=triu(corrcoef(reshapedSpreads))';
corrs(corrs==0)=nan;
mat2Tex(corrs,corrs,{'ISSM','Gibbs','HL','CHL','VoV','LF Ave'},2);

%% Table 2: Zeroing in on the Average anomaly's expected return

fprintf('\n\n\nTable 2 output:\n\n\n');

clear
clc

load organized_results

% must be rows
implist = {'base','opt','optvw'};
samplist = {'insamp','postpub','postpub05'};

coeffs = [];
se = [];
h = {};
for imp = implist
    
    for samp = samplist   

        turn = sumall.(imp{1}).(samp{1}).basic.turnover;
        spreadpaid = sumall.(imp{1}).(samp{1}).basic.spreadpaid;
        gret = sumall.(imp{1}).(samp{1}).g.rbar;
        nret = sumall.(imp{1}).(samp{1}).n.rbar;
        tcost = gret - nret;

        x = 1e4*[gret turn/1e2 spreadpaid tcost nret];

        coeffs = [coeffs; nanmean(x)];
        se = [se; nanstd(x)./sqrt(sum(~isnan(x)))];
        
        switch char(samp)
            case 'insamp'
                h = [h; {'In-Sample'}];
            case 'postpub'
                h = [h; {'Post-Publication'}];
            case 'postpub05'
                h = [h; {'Post-Pub \& Post-2005'}];                
        end

    end
    
end

mat2Tex(coeffs,se,h,0,'(');

%% Table 3: The best expected returns using out-of-sample tests

fprintf('\n\n\nTable 3 output:\n\n\n');

clear
clc

load organized_results

for imp = {'opt','optvw'}

    coeffs = nan(3,1);
    se=nan(3,1);
    h = {};    
    tempstat = sumall.(imp{1});

    for perfi = 1:3

        switch perfi
            case 1
                insampperf = tempstat.insamp.n.rbar;
                perfname = 'nret';
                h=[h; {'Net Return'}];
            case 2
                insampperf = tempstat.insamp.n.SR;
                perfname = 'SR';
                h=[h; {'Net Sharpe'}];
            case 3
                insampperf = -tempstat.insamp.basic.turnover;
                perfname = '-turn';                        
                h=[h; {'1/Turnover'}];
        end

        insampgroup = discretize(insampperf,ptile(insampperf, 0:25:100));

        for groupi = 1:max(insampgroup)
            tempr = 1e4*cropnan(tempstat.postpub05.n.rbar(insampgroup==groupi));    

            coeffs(perfi,groupi+1) = mean(tempr);
            se(perfi,groupi+1) = std(tempr)/sqrt(length(tempr));

        end

    end % perfi

    mat2Tex(coeffs,se,h,1,'(');
    
end % imp

%% Table 4: Empirical Bayes estimates of the best expected returns

fprintf('\n\n\nTable 4 output:\n\n\n');

clear
clc

load organized_results

% Estimation user settings
Nboot = 50;
seed = 1;
plist = [50 70 80 90];

nusrlist = [100 4];

% estimate 
nusr = []; sigsr = []; musr = [];
radj_p50 = []; radj_p60 = []; radj_p70 = []; radj_p80 = []; radj_p90 = [];
impset = [];
vol_radj_p80 = []; vol_radj_p90 = [];
count=1;
for imp = {'opt','optvw'}
    for numui = 1:length(nusrlist)
        
        netstat = sumall.(imp{1}).postpub05.n;        
        tempnusr = nusrlist(numui);
        
        % shrink
        [point,boot] = shrink(netstat,tempnusr,Nboot,plist);

        coeffs(count,:)=[tempnusr sqrt(12)*point.sigsr sqrt(12)*point.musr nan 1e4*point.radj_quant(1:4)];
        se(count,:)=[nan sqrt(12)*std(boot.sigsr) sqrt(12)*std(boot.musr) nan 1e4*std(boot.radj_quant(:,1:4),[],1)];
        
        count=count+1;
    end
end

mat2Tex(coeffs,se,{},2,'(');
mat2Tex(coeffs,se,{},1,'(');


%% Table 5: Combination strategies

fprintf('\n\n\nTable 5 output:\n\n\n');

clear
clc

load res_combo_strats

nstrats=size(resISew,2);
ncombos=size(resISew,1);     


xretIS=nan(ncombos,nstrats);
xretOS=nan(ncombos,nstrats);

txretIS=nan(ncombos,nstrats);
txretOS=nan(ncombos,nstrats);

netxretIS=nan(ncombos,nstrats);
netxretOS=nan(ncombos,nstrats);

tnetxretIS=nan(ncombos,nstrats);
tnetxretOS=nan(ncombos,nstrats);

toIS=nan(ncombos,nstrats);
toOS=nan(ncombos,nstrats);


for i=1:size(resOSew,1)
    for j=1:size(resOSew,2)
        xretIS(i,j)=resISew(i,j).xret(end); 
        txretIS(i,j+nstrats)=resISvw(i,j).xret(end);
        xretOS(i,j)=resOSew(i,j).xret(end); 
        xretOS(i,j+nstrats)=resOSvw(i,j).xret(end);        

        txretIS(i,j)=resISew(i,j).txret(end); 
        txretIS(i,j+nstrats)=resISvw(i,j).txret(end);
        txretOS(i,j)=resOSew(i,j).txret(end); 
        txretOS(i,j+nstrats)=resOSvw(i,j).txret(end);        

        netxretIS(i,j)=resISew(i,j).netxret; 
        netxretIS(i,j+nstrats)=resISvw(i,j).netxret;
        netxretOS(i,j)=resOSew(i,j).netxret; 
        netxretOS(i,j+nstrats)=resOSvw(i,j).netxret;        

        tnetxretIS(i,j)=resISew(i,j).tnetxret; 
        tnetxretIS(i,j+nstrats)=resISvw(i,j).tnetxret;
        tnetxretOS(i,j)=resOSew(i,j).tnetxret; 
        tnetxretOS(i,j+nstrats)=resOSvw(i,j).tnetxret;        

        toIS(i,j)=resISew(i,j).turnover; 
        toIS(i,j+nstrats)=resISvw(i,j).turnover;
        toOS(i,j)=resOSew(i,j).turnover; 
        toOS(i,j+nstrats)=resOSvw(i,j).turnover;        
        
    end
end

a=[];
tA=[];

for r=1:4
    c=find(netxretIS(r,:)==max(netxretIS(r,:)));
    
    a=[a; xretIS(r,1) toIS(r,1) netxretIS(r,1)];
    tA=[tA; txretIS(r,1) nan tnetxretIS(r,1)];
    
    a=[a; xretOS(r,1) toOS(r,1) netxretOS(r,1)];
    tA=[tA; txretOS(r,1) nan tnetxretOS(r,1)];

    a=[a; xretIS(r,c) toIS(r,c) netxretIS(r,c)];
    tA=[tA; txretIS(r,c) nan tnetxretIS(r,c)];
    
    a=[a; xretOS(r,c) toOS(r,c) netxretOS(r,c)];
    tA=[tA; txretOS(r,c) nan tnetxretOS(r,c)];
end

a=100*[a(1:8,:) nan(8,1) a(9:end,:)];
tA=[tA(1:8,:) nan(8,1) tA(9:end,:)];

h={'1985-2005','2006-2020','1985-2005','2006-2020','1985-2005','2006-2020','1985-2005','2006-2020'};

mat2Tex(a,a./tA,h,0,'bracketsType','(');


%% Table 6: Reconciliation with NMV

fprintf('\n\n\nTable 6 output:\n\n\n');

clear
clc

panelA=[];
panelB=[];

load res_nmv_reconciliation1

for i=1:length(res)
    res(i).gross_ret=res(i).xret(end);
    res(i).net_ret=res(i).netxret(end);
end

panelA=[panelA [mean([res(find([res.turnover]<0.10)).gross_ret]) ...
    mean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).gross_ret]) ...
    mean([res(find([res.turnover]'>=0.50)).gross_ret]) ...
    mean([res.gross_ret])]'];

panelB=[panelB [mean([res(find([res.turnover]'<0.10)).netxret]) ...
    mean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).netxret]) ...
    mean([res(find([res.turnover]'>=0.50)).netxret]) ...
    mean([res.netxret])]'];


res(1).postPubPost2005Xret=nan;
res(1).postPubPost2005Netret=nan;

for i=1:length(res)
    postPubPost2005Start=find(res(i,1).dates>100*res(i,1).pubYear+12 & res(i,1).dates>200512,1,'first');
    res(i).postPubPost2005Xret=100*nanmean(res(i).pret(postPubPost2005Start:end,end));
    res(i).postPubPost2005Netret=100*nanmean(res(i).netpret(postPubPost2005Start:end));    
end

% Post pub & post 2005
panelA=[panelA [mean([res(find([res.turnover]'<0.10)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Xret]) ...
    nanmean([res.postPubPost2005Xret])]'];

panelB=[panelB [mean([res(find([res.turnover]'<0.10)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Netret]) ...
    nanmean([res.postPubPost2005Netret])]'];


load res_nmv_reconciliation2
mitigatedIndicator=~(strcmp(signal_names,'ChangeInRecommendation') | strcmp(signal_names,'NumEarnIncrease') | strcmp(signal_names,'RDcap') | strcmp(anomaly_summary.Cat_Form,'discrete'));
res=res(mitigatedIndicator,:);
anomaly_summary=anomaly_summary(mitigatedIndicator,:);
signal_names=signal_names(mitigatedIndicator,:);


for i=1:length(res)
    res(i).gross_ret=res(i).xret(end);
    res(i).net_ret=res(i).netxret(end);
end

panelA=[panelA [mean([res(find([res.turnover]<0.10)).gross_ret]) ...
    mean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).gross_ret]) ...
    mean([res(find([res.turnover]'>=0.50)).gross_ret]) ...
    mean([res.gross_ret])]'];

panelB=[panelB [mean([res(find([res.turnover]'<0.10)).netxret]) ...
    mean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).netxret]) ...
    mean([res(find([res.turnover]'>=0.50)).netxret]) ...
    mean([res.netxret])]'];



load res_nmv_reconciliation3
stayIndicator=~strcmp(signal_names,'IO_ShortInterest');
res=res(stayIndicator,:);
anomaly_summary=anomaly_summary(stayIndicator,:);
signal_names=signal_names(stayIndicator,:);
res(1).postPubPost2005Xret=nan;
res(1).postPubPost2005Netret=nan;

for i=1:length(res)
    postPubPost2005Start=find(res(i,1).dates>100*anomaly_summary.Year(i)+12 & res(i,1).dates>200512,1,'first');
    res(i).postPubPost2005Xret=100*nanmean(res(i).pret(postPubPost2005Start:end,end));
    res(i).postPubPost2005Netret=100*nanmean(res(i).netpret(postPubPost2005Start:end));    
end

% Post pub & post 2005
panelA=[panelA [nanmean([res(find([res.turnover]'<0.10)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Xret]) ...
    nanmean([res.postPubPost2005Xret])]'];

panelB=[panelB [nanmean([res(find([res.turnover]'<0.10)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Netret]) ...
    nanmean([res.postPubPost2005Netret])]'];


load res_unmitigated
res=res(:,3);
stayIndicator=~strcmp(signal_names,'IO_ShortInterest');
res=res(stayIndicator,:);
anomaly_summary=anomaly_summary(stayIndicator,:);
signal_names=signal_names(stayIndicator,:);


res(1).postPubPost2005Xret=nan;
res(1).postPubPost2005Netret=nan;

for i=1:length(res)
    postPubPost2005Start=find(res(i,1).dates>=100*anomaly_summary.Year(i)+12 & res(i,1).dates>200512,1,'first');
    res(i).postPubPost2005Xret=100*nanmean(res(i).pret(postPubPost2005Start:end,end));
    res(i).postPubPost2005Netret=100*nanmean(res(i).netpret(postPubPost2005Start:end));    
end

% Post pub & post 2005
panelA=[panelA [nanmean([res(find([res.turnover]'<0.10)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Xret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Xret]) ...
    nanmean([res.postPubPost2005Xret])]'];

panelB=[panelB [nanmean([res(find([res.turnover]'<0.10)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.10 & [res.turnover]'<0.50)).postPubPost2005Netret]) ...
    nanmean([res(find([res.turnover]'>=0.50)).postPubPost2005Netret]) ...
    nanmean([res.postPubPost2005Netret])]'];


columnNames=array2table(reshape({'Anomaly','NV (2016)','NV (2016)','NV (2016)','CZ (2020)','CZ (2020)',...
            'Strategy Construction','VW NYSE','VW NYSE','VW NYSE','Orig','Orig',...
            'Tcosts','Gibbs','Gibbs','Gibbs','Gibbs','Combined', ...
            'Sample','Full','Post-pub+2005','Full','Post-pub+2005','Post-pub+2005'},6,4)','VariableNames',{'Var1','Col1','Col2','Col3','Col4','Col5'});
disp(columnNames);



h={'Low','Medium','High','All'}';
mat2Tex(100*panelA,100*panelA,h,0);
mat2Tex(100*panelB,100*panelB,h,0);


%% Table 7: Reconcilliation with institutional trading cost studies 

fprintf('\n\n\nTable 7 output:\n\n\n');

clear
clc

load organized_results

anomilist = [ ...
    find(strcmp('Size',anomaly_summary.Acronym)) ...
    find(strcmp('BMdec',anomaly_summary.Acronym)) ...
    find(strcmp('Mom12m',anomaly_summary.Acronym)) ...
    ];

coeffs=[];
se=[];
for sampi = 1:2
    for anomi = anomilist
        switch sampi
            case 1
                t = (resall.opt(anomi).dates >= 200601);
            case 2
                t = (resall.opt(anomi).dates >= 199801) ...
                    & (resall.opt(anomi).dates <= 201312);
        end

        gret = 1e4*resall.opt(anomi).pret(t,end);
        nret = 1e4*resall.opt(anomi).netpret(t);
        
        coeffs=[coeffs; mean(gret); mean(nret)];
        se=[se; std(gret)/sqrt(sum(t)); std(nret)/sqrt(sum(t))];
    end       
end

coeffs=reshape(coeffs,6,2);
se=reshape(se,6,2);

% Add the FIM, BLNR, and PW results
rightHandSidePanelCoeffs=[66.5 44.8 16.3; 
                    54.3 43.4 18.3;
                    40.5 25.3 45.3;
                    29.3 23.3 19.3;
                    18.8 45.9 50.1;
                    -6.4 23.1 14.4];
rightHandSidePanelSE=[22.1 28.1 20.2;                
                    21.9 nan 19.9;
                    36.2 31.2 23.4;
                    36.6 nan 23.2;
                    47.1 47.4 31.3;
                    45.8 nan 32.0];
coeffs=[coeffs nan(6,1) rightHandSidePanelCoeffs];
se=[se nan(6,1) rightHandSidePanelSE];

h={'Gross','Net','Gross','Net','Gross','Net'};
mat2Tex(coeffs,se,h,1,'(');

%% Timekeeping

diary off

fprintf('\n\n\n\nDone with table printing @ %s\n\n\n',char(datetime('now')));
