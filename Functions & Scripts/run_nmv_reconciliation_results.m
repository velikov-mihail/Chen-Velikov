% Novy-Marx and Velikov (2016) replication

clear
clc

load dates
load me
load ret
load NYSE
load gibbs_filled

fprintf('\n NMV reconcilitaion results part 1. Run started at: ');
disp(datetime('now'));
tic;

[anoms23, labels23]=getAnomalySignals('novyMarxVelikovAnomalies.csv',1,2);

startDate=[1 1 1 1 1 1 1 1 2 2 2 1 1 1 1 2 2 1 1 1 1 1 1]; % 1 - 196306, 2 - 197306
startDate(startDate==1)=196306;
startDate(startDate==2)=197306;

for i=1:size(anoms23,3)
    ind=makeUnivSortInd(anoms23(:,:,i),10,NYSE);
    res(i,1) = runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'timePeriod',[startDate(i) 202012],'plotFigure',0,'printResults',0); 
end

pubYear=[1993 2013 1993 2014 1996 2008 2008 2000 2008 2010 2008 2014 2014 2006 1993 1984 2008 1999 2014 2016 1993 2011 2016]; % Loosely ad-hoc determined publication years
for i=1:size(anoms23,3)
    res(i,1).label=labels23(i);
    res(i,1).pubYear=pubYear(i);
end

save Results\res_nmv_reconciliation1 res

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;

%% CZ anomalies, NMV strategy constrution, gibbs tcosts

clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load gibbs_filled

fprintf('\n Now working on NMV reconcilitaion results part 2. Run started at: ');
disp(datetime('now'));
tic;

% Get the anomalies; this step requires a lot of memory
[anoms,labels,anomaly_summary]=getChenZimmermanAnomalies();

JuneIndicator=(dates-100*floor(dates/100)==6);
JuneDecemberIndicator=(dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator=(dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

signal_names=anomaly_summary.Acronym;

for i=1:height(anomaly_summary)
    
    % Assign the current anomaly signal to a matrix called 'var'
    r=find(strcmp(labels,anomaly_summary.Acronym(i)));    
    var=anoms(:,:,r);
        
    % Figure out the rebalancing period
    if anomaly_summary.PortfolioPeriod(i)==12 | anomaly_summary.PortfolioPeriod(i)==36
        var(~JuneIndicator,:)=nan;
    elseif anomaly_summary.PortfolioPeriod(i)==6
        var(~JuneDecemberIndicator,:)=nan;
    elseif anomaly_summary.PortfolioPeriod(i)==3
        var(~quarterEndIndicator,:)=nan;
    end
        
    % Check for sample filters
    if ~strcmp(anomaly_summary.Filter(i),'')
        filterString=regexprep(anomaly_summary.Filter(i),' ','');
        if contains(filterString,'abs(prc)>1')
            var(abs(prc)<=1)=nan;  
        end
        if contains(filterString,'abs(prc)>5')
            var(abs(prc)<=5)=nan;  
        end
        if contains(filterString,'exchcd%in%c(1,2)')
            var(exchcd>2)=nan;
        end
        if contains(filterString,'exchcd==1')
            var(exchcd>1)=nan;
        end
        if contains(filterString,'me>me_nyse20')
            indME5=makeUnivSortInd(me,5,NYSE);
            var(indME5==1)=nan;
        end
    end
                        
    % Determine the sample
    if isnan(anomaly_summary.StartMonth(i))
        sampleStart=anomaly_summary.SampleStartYear(i)*100+6;    
    else
        sampleStart=anomaly_summary.SampleStartYear(i)*100+anomaly_summary.StartMonth(i);
    end
    sampleEnd=dates(min(find(sum(isfinite(var),2)>0,1,'last')+anomaly_summary.PortfolioPeriod(i),length(dates)));
    dts=[sampleStart sampleEnd];
    
    % Determine the sorting breakpoints  
    if strcmp(anomaly_summary.Acronym(i),'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i),'NumEarnIncrease')
        bpts = [];
        for j = 1:4
            bpts = [bpts j*100/5];
        end
        bpts=prctile(var,bpts,2);        
        ind=1*(var<=repmat(bpts(:,1),1,size(var,2)))+2*(var>=repmat(bpts(:,end),1,size(var,2)));        
    elseif strcmp(anomaly_summary.Acronym(i),'RDcap')
        ind=1*(var==0)+2*(var>0);
    elseif strcmp(anomaly_summary.Cat_Form(i),'discrete')       
        uVals=unique(var);
        uVals(isnan(uVals))=[];
        ind=1*(var==min(uVals))+2*(var==max(uVals));
    else
        ind = makeUnivSortInd(var,10,NYSE);  
    end
    res(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,size(anoms,3));
end

save Results\res_nmv_reconciliation2 res anomaly_summary signal_names

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;

%% CZ anomalies, original strategy constrution, gibbs tcosts

clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load gibbs_filled

fprintf('\n NMV reconcilitaion results part 3. Run started at: ');
disp(datetime('now'));
tic;

% Get the anomalies; this step requires a lot of memory
[anoms,labels,anomaly_summary]=getChenZimmermanAnomalies();

JuneIndicator=(dates-100*floor(dates/100)==6);
JuneDecemberIndicator=(dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator=(dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

signal_names=anomaly_summary.Acronym;

for i=1:height(anomaly_summary)
    
    % Assign the current anomaly signal to a matrix called 'var'
    r=find(strcmp(labels,anomaly_summary.Acronym(i)));    
    var=anoms(:,:,r);
        
    % Figure out the rebalancing period
    if anomaly_summary.PortfolioPeriod(i)==12 || anomaly_summary.PortfolioPeriod(i)==36
        var(~JuneIndicator,:)=nan;
    elseif anomaly_summary.PortfolioPeriod(i)==6
        var(~JuneDecemberIndicator,:)=nan;
    elseif anomaly_summary.PortfolioPeriod(i)==3
        var(~quarterEndIndicator,:)=nan;
    end
        
    % Check for sample filters
    if ~strcmp(anomaly_summary.Filter(i),'')
        filterString=regexprep(anomaly_summary.Filter(i),' ','');
        if contains(filterString,'abs(prc)>1')
            var(abs(prc)<=1)=nan;  
        end
        if contains(filterString,'abs(prc)>5')
            var(abs(prc)<=5)=nan;  
        end
        if contains(filterString,'exchcd%in%c(1,2)')
            var(exchcd>2)=nan;
        end
        if contains(filterString,'exchcd==1')
            var(exchcd>1)=nan;
        end
        if contains(filterString,'me>me_nyse20')
            indME5=makeUnivSortInd(me,5,NYSE);
            var(indME5==1)=nan;
        end
    end
                        
    % Determine the sample
    if isnan(anomaly_summary.StartMonth(i))
        sampleStart=anomaly_summary.SampleStartYear(i)*100+6;    
    else
        sampleStart=anomaly_summary.SampleStartYear(i)*100+anomaly_summary.StartMonth(i);
    end
    sampleEnd=dates(min(find(sum(isfinite(var),2)>0,1,'last')+anomaly_summary.PortfolioPeriod(i),length(dates)));
    dts=[sampleStart sampleEnd];    
    
    % Determine the sorting breakpoints 
    if strcmp(anomaly_summary.Acronym(i),'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i),'NumEarnIncrease')
        bpts = [];
        for j = 1:4
            bpts = [bpts j*100/5];
        end
        bpts=prctile(var,bpts,2);        
        ind=1*(var<=repmat(bpts(:,1),1,size(var,2)))+2*(var>=repmat(bpts(:,end),1,size(var,2)));        
    elseif strcmp(anomaly_summary.Acronym(i),'RDcap')
        ind=1*(var==0)+2*(var>0);
    elseif strcmp(anomaly_summary.Cat_Form(i),'discrete')       
        uVals=unique(var);
        uVals(isnan(uVals))=[];
        ind=1*(var==min(uVals))+2*(var==max(uVals));
    else
        if isnan(anomaly_summary.LSQuantile(i))
            nptfs=5;
        else
            nptfs=round(1/anomaly_summary.LSQuantile(i));
        end
        ind = makeUnivSortInd(var,nptfs);  
        if strcmp(anomaly_summary.QuantileFilter(i),'NYSE')
            ind = makeUnivSortInd(var,nptfs,NYSE);  
        end
    end
    if strcmp(anomaly_summary.StockWeight(i),'VW')
        res(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 
    else
        res(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);
    end
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,size(anoms,3));
end

save Results\res_nmv_reconciliation3 res anomaly_summary signal_names

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;
