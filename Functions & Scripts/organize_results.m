% using matlab 2021a

clear
clc

% ==== IMPORT DATA ====
raw_bh = load('res_buy_hold_all.mat');
raw_lc = load('res_mitigated_low_cost.mat');
raw_reb = load('res_reduced_reb_freq.mat');
raw_unmit = load('res_unmitigated.mat');

anomaly_summary = raw_unmit.anomaly_summary;

% Clean up anomaly characteristics data
anomaly_summary.PortfolioPeriod(isnan(anomaly_summary.PortfolioPeriod))=1;
anomaly_summary.PortfolioPeriod(anomaly_summary.PortfolioPeriod>12)=12;

anomaly_summary.StockWeight(strcmp(anomaly_summary.StockWeight,'NA'))= {'EW'};
anomaly_summary.StockWeight=lower(anomaly_summary.StockWeight);

% res_all is 205 signals by N implementations
res_too_many = [
    raw_bh.res_ew(:,2:end) ...
    raw_bh.res_vw(:,2:end) ...
    raw_lc.res_ew ...
    raw_lc.res_vw ...    
    raw_reb.res_ew ...
    raw_reb.res_vw ...        
    ];

% remove IO_ShortInterest
% it has too few stocks, only a handful per month in each leg in the in-sample period
ibad = find(strcmp(anomaly_summary.Acronym,'IO_ShortInterest'));
res_too_many(ibad,:) = [];
anomaly_summary(ibad,:) = [];
raw_unmit.anomaly_summary(ibad,:)=[];
raw_unmit.res(ibad,:)=[];
raw_unmit.signal_names(ibad,:)=[];


implist = {
    'ew' 'bh2025'
    'ew' 'bh2030'
    'ew' 'bh2035'
    'ew' 'bh2040'    
    'ew' 'bh2045'        
    'ew' 'bh2050'
    'vw' 'bh1020'
    'vw' 'bh1030'
    'vw' 'bh1040'
    'vw' 'bh1050'
    'ew' 'lc1'
    'ew' 'lc2'
    'vw' 'lc1'
    'vw' 'lc2'    
    'ew' 'reb01'
    'ew' 'reb03'
    'ew' 'reb06'
    'ew' 'reb12'
    'vw' 'reb01'
    'vw' 'reb03'
    'vw' 'reb06'
    'vw' 'reb12'}';


% make a list of rebalancing periods (not sure why this is so hard)
f = @(s) strcmp(s(1:3),'reb');
f2 = @(s) str2num(s(end-1:end));
ireb = cellfun(f, implist(2,:));
reblist = nan(size(implist(2,:)));
reblist(ireb) = cell2mat(cellfun(f2, implist(2,ireb), 'UniformOutput', 0));


% create resall
% resall consists of organized 205  res structs

clear resall
for anomi = 1:size(res_too_many,1)    
    
    % create months_since_pub for adding on later    
    dates = res_too_many(anomi,1).dates;    
    year = floor(dates/100);
    month = dates-year*100;    
    tcal = (year)*12 + month;    
    months_since_pub = tcal - anomaly_summary.Year(anomi)*12-12;
    
    % find mean ret in samp for optimizing later
    obj = nan(1,size(res_too_many,2));    
    tok = res_too_many(anomi,1).dates >= anomaly_summary.SampleStartYear(anomi)*100+anomaly_summary.StartMonth(anomi) ...
            & res_too_many(anomi,1).dates <= anomaly_summary.SampleEndYear(anomi)*100+12;
    
    for impi = 1:size(res_too_many,2)                               
        obj(1,impi) = nanmean(res_too_many(anomi,impi).netpret(tok));
    end % if continuous
        
    temp = raw_unmit.res(anomi,3);
    temp.sweight = char(anomaly_summary.StockWeight(anomi));
    temp.imp     = anomaly_summary.PortfolioPeriod(anomi);
    temp.months_since_pub = months_since_pub;
    resall.base(anomi,1) = temp;                        
    
    % find imp indexes for reb that is not mitigation (and skip)
    impskip = reblist <= anomaly_summary.PortfolioPeriod(anomi);       

    % save optimized over all implementations
    tempobj = obj;  
    tempobj(impskip) = -inf;    
    [~,impstar] = max(tempobj);        
    temp = res_too_many(anomi,impstar);
    temp.sweight = implist{1,impstar};
    temp.imp     = reblist(1,impstar);    
    temp.months_since_pub = months_since_pub;
    resall.opt(anomi,1) = temp;    
    
    % save optimized over vw implementations
    tempobj = obj;
    tempobj(~strcmp(implist(1,:),'vw')) = -inf;
    tempobj(impskip) = -inf;
    [~,impstar] = max(tempobj);    
    temp = res_too_many(anomi,impstar);
    temp.sweight = implist{1,impstar};
    temp.imp     = reblist(1,impstar);        
    temp.months_since_pub = months_since_pub;
    resall.optvw(anomi,1) = temp;        
    
end % for res entry i


% SUMMARY STATS summary stats (sumall)

impcatlist = fieldnames(resall)';


clear sumall

% first create sumall.set a (nstrat x 1) table of implmenetations / settings

for impcat = impcatlist
    clear signalname sweight imp 
    for anomi = 1:length(resall.(impcat{1}))
        tempts = resall.(impcat{1})(anomi);

        signalname{anomi,1} = anomaly_summary.Acronym(anomi);
        sweight{anomi,1} = tempts.sweight;
        imp{anomi,1} = tempts.imp;    

        sampstart(anomi,1) = anomaly_summary.SampleStartYear(anomi)*100+anomaly_summary.StartMonth(anomi);
        sampend(anomi,1) = anomaly_summary.SampleEndYear(anomi)*100+12;    
        pubdate(anomi,1) = anomaly_summary.Year(anomi)*100+12;            
    end    
    
    sumall.(impcat{1}).set = table(signalname, sweight, imp, sampstart, sampend, pubdate);
end
    

% add performance stats

for impcat = impcatlist
    
    for samptype = {'insamp','postpub','postpub05' ...
            ,'from85to14','from85to05','from06to14','from06to20','from10to14'}
        for rettype = {'g','n'}
            
            [rbar,vol,T, turnover,spreadpaid] = nanall(size(anomaly_summary,1),1);
            
            for anomi = 1:size(anomaly_summary,1)
                
                tempts = resall.(impcat{1})(anomi);               

                % find dates
                switch samptype{1}
                    case 'insamp'
                        sampstart = sumall.(impcat{1}).set(anomi,:).sampstart;
                        sampend = sumall.(impcat{1}).set(anomi,:).sampend;
                    case 'postpub'
                        sampstart = sumall.(impcat{1}).set(anomi,:).pubdate;
                        sampend = inf;
                    case 'postpub05'
                        sampstart = max(sumall.(impcat{1}).set(anomi,:).pubdate, 200601);
                        sampend = inf;
                    case 'from85to14'
                        sampstart = 198501; sampend = 201412;
                    case 'from85to05'
                        sampstart = 198501; sampend = 200512;
                    case 'from06to14'
                        sampstart = 200601; sampend   = 201412;                        
                    case 'from06to20'
                        sampstart = 200601; sampend   = 202012;                                                
                    case 'from10to14'
                        sampstart = 201001; sampend   = 201412;
                end

                t = tempts.dates >= sampstart & tempts.dates <= sampend;                                        
                
                % choose ret type
                switch rettype{1}
                    case 'g'
                        ret = tempts.pret(t,end);
                    case 'n'
                        ret = tempts.netpret(t,end);
                end
                
                % find stats
                rbar(anomi) = nanmean(ret);
                vol(anomi) = nanstd(ret);
                T(anomi) = sum(~isnan(ret));
                
                % turnover should be two sided (max 200% for LS)
                turnover(anomi) = nanmean(nanmean(tempts.toTS(t,:),2));
                spreadpaid(anomi) = nanmean(...
                    (tempts.pret(t,end) - tempts.netpret(t,end)) ...
                    ./nanmean(tempts.toTS(t,:),2)...
                    );
                
            end
            
            tstat = rbar./vol.*(T.^(0.5));
            SR = rbar./vol;
            sumall.(impcat{1}).(samptype{1}).(rettype{1}) = table(signalname,rbar,vol,tstat,SR);
            
        end % for rettype
        
        sumall.(impcat{1}).(samptype{1}).basic = table(signalname,T,turnover,spreadpaid);
        
    end % for samptype
end % for impcat

% save cleandat for future use
save('Results/organized_results','resall','sumall','anomaly_summary')
