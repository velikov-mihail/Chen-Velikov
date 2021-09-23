
clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load tcosts_cv

fprintf('\n Now working on unmitigated results. Run started at: ');
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
                        
    % Determine the sample period
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
    
    % Run equal- and value-weighted sorts
    res(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);      
    res(i,2)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0);  
    % The third column in the structure stores the weighting from the
    % original publication
    if strcmp(anomaly_summary.StockWeight(i),'VW')
        res(i,3)=res(i,2);
    else
        res(i,3)=res(i,1);
    end
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,size(anoms,3));
end

% Store the structure with the results as well as the summary tables
save Results\res_unmitigated res anomaly_summary signal_names

% Time-keeping
fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;
