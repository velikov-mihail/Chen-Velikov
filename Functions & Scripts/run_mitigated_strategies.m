%% Does all 5 VW and 7 EW buy/hold strategies

clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load tcosts_cv

fprintf('\n Now working on mitigated buy/hold results. Run started at: ');
disp(datetime('now'));
tic;

% Get the anomalies; this step requires a lot of memory
[anoms,labels,anomaly_summary]=getChenZimmermanAnomalies();

JuneIndicator=(dates-100*floor(dates/100)==6);
JuneDecemberIndicator=(dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator=(dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

% We are not doing buy/hold for the discrete signals plus a few other problematic ones
mitigatedIndicator=~(strcmp(anomaly_summary.Acronym,'ChangeInRecommendation') | strcmp(anomaly_summary.Acronym,'NumEarnIncrease') | strcmp(anomaly_summary.Acronym,'RDcap') | strcmp(anomaly_summary.Cat_Form,'discrete'));

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
    
    if mitigatedIndicator(i)     
        
        % Equal-weighted name breaks (20/20)
        ind = makeUnivSortInd(var,5);    
        ind=1*(ind==1)+2*(ind==5);
        res_ew(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);    

        % Equal-weighted name breaks (20/25)
        ind = makeUnivSortInd(var,20);    
        ind(ind==2 | ind==3 | ind==4)=1;
        ind(ind==19 | ind==18| ind==17)=20;
        ind=make_sS(ind,5);
        ind=1*(ind==1)+2*(ind==20);
        res_ew(i,2)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);    
        
        
        % Equal-weighted name breaks (20/30)
        ind = makeUnivSortInd(var,10);    
        ind(ind==2)=1;
        ind(ind==9)=10;
        ind=make_sS(ind,3);
        ind=1*(ind==1)+2*(ind==10);
        res_ew(i,3)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);   

        % Equal-weighted name breaks (20/35)
        ind = makeUnivSortInd(var,20);    
        ind(ind==2 | ind==3 | ind==4)=1;
        ind(ind==19 | ind==18| ind==17)=20;
        ind=make_sS(ind,7);
        ind=1*(ind==1)+2*(ind==20);
        res_ew(i,4)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);           
        
        % Equal-weighted name breaks (20/40)
        ind = makeUnivSortInd(var,10);    
        ind(ind==2)=1;
        ind(ind==9)=10;
        ind=make_sS(ind,4);
        ind=1*(ind==1)+2*(ind==10);
        res_ew(i,5)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);    

        % Equal-weighted name breaks (20/45)
        ind = makeUnivSortInd(var,20);    
        ind(ind==2 | ind==3 | ind==4)=1;
        ind(ind==19 | ind==18| ind==17)=20;
        ind=make_sS(ind,9);
        ind=1*(ind==1)+2*(ind==20);
        res_ew(i,6)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);   
        
        % Equal-weighted name breaks (20/50)
        ind = makeUnivSortInd(var,10);    
        ind(ind==2)=1;
        ind(ind==9)=10;
        ind=make_sS(ind,5);
        ind=1*(ind==1)+2*(ind==10);
        res_ew(i,7)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);    
        
        % Value-weighted NYSE breaks (10/10)
        ind = makeUnivSortInd(var,10,NYSE);    
        ind=1*(ind==1)+2*(ind==10);
        res_vw(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0);
    
         % Value-weighted NYSE breaks (10/20)
        ind = makeUnivSortInd(var,10,NYSE);    
        ind=make_sS(ind,2);
        ind=1*(ind==1)+2*(ind==10);
        res_vw(i,2)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 
   
         % Value-weighted NYSE breaks (10/30)
        ind = makeUnivSortInd(var,10,NYSE);    
        ind=make_sS(ind,3);
        ind=1*(ind==1)+2*(ind==10);
        res_vw(i,3)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 

         % Value-weighted NYSE breaks (10/40)
        ind = makeUnivSortInd(var,10,NYSE);    
        ind=make_sS(ind,4);
        ind=1*(ind==1)+2*(ind==10);
        res_vw(i,4)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 
    
    
         % Value-weighted NYSE breaks (10/50)
        ind = makeUnivSortInd(var,10,NYSE);    
        ind=make_sS(ind,5);
        ind=1*(ind==1)+2*(ind==10);
        res_vw(i,5)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0); 
   
    end
    
    clear var
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,length(labels));
end


% Replace the ones we are not mitigating with the unmitigated results
load res_unmitigated
for i=1:length(mitigatedIndicator)
    if ~mitigatedIndicator(i)
        res_ew(i,:)=repmat(res(i,1),1,size(res_ew,2));
        res_vw(i,:)=repmat(res(i,2),1,size(res_vw,2));
    end
end

% Store the structure with the results as well as the summary tables
save Results\res_buy_hold_all res_ew res_vw anomaly_summary signal_names

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;

%% Does all 4 VW and 4 EW reduced rebalancing frequency strategies

clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load tcosts_cv

fprintf('\n Now working on mitigated reduced rebalancing results. Run started at: ');
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
            
    % Loop through monthly, quarterly, semi-annual, and annual rebalancing
    for j=1:4
        if j==1
            temp_var=var;
        elseif j==2
            temp_var=var;
            temp_var(~quarterEndIndicator,:)=nan;
        elseif j==3
            temp_var=var;
            temp_var(~JuneDecemberIndicator,:)=nan;            
        else                    
            temp_var=var;
            temp_var(~JuneIndicator,:)=nan;                    
        end
        
       % Determine the sorting breakpoints 
       if strcmp(anomaly_summary.Acronym(i),'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i),'NumEarnIncrease')
            bpts = [];
            for j = 1:4
                bpts = [bpts j*100/5];
            end
            bpts=prctile(temp_var,bpts,2);        
            ind=1*(temp_var<=repmat(bpts(:,1),1,size(temp_var,2)))+2*(temp_var>=repmat(bpts(:,end),1,size(temp_var,2)));        
        elseif strcmp(anomaly_summary.Acronym(i),'RDcap')
            ind=1*(temp_var==0)+2*(temp_var>0);
        elseif strcmp(anomaly_summary.Cat_Form(i),'discrete')       
            uVals=unique(temp_var);
            uVals(isnan(uVals))=[];
            ind=1*(temp_var==min(uVals))+2*(temp_var==max(uVals));
        else
            if isnan(anomaly_summary.LSQuantile(i))
                nptfs=5;
            else
                nptfs=round(1/anomaly_summary.LSQuantile(i));
            end
            ind = makeUnivSortInd(temp_var,nptfs);  
            if strcmp(anomaly_summary.QuantileFilter(i),'NYSE')
                ind = makeUnivSortInd(temp_var,nptfs,NYSE);  
            end
        end
        res_ew(i,j)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);      
        res_vw(i,j)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0);
    end
    
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,size(anoms,3));
end

% Replace the ones where we had issues with the mitigation with the unmitigated results
load res_unmitigated
for i=1:size(res_ew,1)
    for j=1:size(res_ew,2)
        if isempty(res_ew(i,j).xret)
            res_ew(i,j)=res(i,1);
        end
        if isempty(res_vw(i,j).xret)
            res_vw(i,j)=res(i,2);
        end
    end
end

save Results\res_reduced_reb_freq res_ew res_vw anomaly_summary signal_names

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;


%% Test low cost strategies

clear
clc

load prc
load ret
load me
load dates
load exchcd
load NYSE
load tcosts_cv

fprintf('\n Now working on mitigated low cost universe results. Run started at: ');
disp(datetime('now'));
tic;

% Get the anomalies; this step requires a lot of memory
[anoms,labels,anomaly_summary]=getChenZimmermanAnomalies();

JuneIndicator=(dates-100*floor(dates/100)==6);
JuneDecemberIndicator=(dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator=(dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

% Create the low cost indicators 
indME10=makeUnivSortInd(me,10,NYSE);
indLowCostTertile=nan(size(ret));
indLowCostHalf=nan(size(ret));

% Low cost tertile/half within each NYSE market cap decile
for i=1:10
    temp=tcosts;
    temp(indME10~=i)=nan;
    ind=makeUnivSortInd(temp,3);
    indLowCostTertile(ind==1)=1;
    ind=makeUnivSortInd(temp,2);
    indLowCostHalf(ind==1)=1;
end


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
    
    % Low cost tertile
    var1=var.*indLowCostTertile; 
    
    if strcmp(anomaly_summary.Acronym(i),'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i),'NumEarnIncrease')
        bpts = [];
        for j = 1:4
            bpts = [bpts j*100/5];
        end
        bpts=prctile(var1,bpts,2);        
        ind=1*(var1<=repmat(bpts(:,1),1,size(var,2)))+2*(var1>=repmat(bpts(:,end),1,size(var,2)));        
    elseif strcmp(anomaly_summary.Acronym(i),'RDcap')
        ind=1*(var1==0)+2*(var1>0);
    elseif strcmp(anomaly_summary.Cat_Form(i),'discrete')       
        uVals=unique(var);
        uVals(isnan(uVals))=[];
        ind=1*(var==min(uVals))+2*(var==max(uVals));
    else
        ind=makeUnivSortInd(var1,5);  
    end
    res_ew(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);      
    res_vw(i,1)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0);      
    
    % Low cost half
    var2=var.*indLowCostHalf; 
    
    if strcmp(anomaly_summary.Acronym(i),'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i),'NumEarnIncrease')
        bpts = [];
        for j = 1:4
            bpts = [bpts j*100/5];
        end
        bpts=prctile(var,bpts,2);        
        ind=1*(var<=repmat(bpts(:,1),1,size(var,2)))+2*(var>=repmat(bpts(:,end),1,size(var,2)));        
    elseif strcmp(anomaly_summary.Acronym(i),'RDcap')
        ind=1*(var1==0)+2*(var1>0);
    elseif strcmp(anomaly_summary.Cat_Form(i),'discrete')       
        uVals=unique(var);
        uVals(isnan(uVals))=[];
        ind=1*(var==min(uVals))+2*(var==max(uVals));
    else
        ind=makeUnivSortInd(var1,5);  
    end
    res_ew(i,2)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',dts,'plotFigure',0,'printResults',0);      
    res_vw(i,2)=runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',dts,'plotFigure',0,'printResults',0);      
    
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n',char(labels(r)),i,size(anoms,3));
end

% Replace the ones where we had issues with the mitigation with the unmitigated results
load res_unmitigated
for i=1:size(res_ew,1)
    for j=1:size(res_ew,2)
        if isempty(res_ew(i,j).xret) || isempty(res_vw(i,j).xret)
            res_ew(i,j)=res(i,1);
            res_vw(i,j)=res(i,2);
        end
    end
end

save Results\res_mitigated_low_cost res_ew res_vw anomaly_summary signal_names

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;
