
clear
clc

load ret
load dates
load nyse
load me
load tcosts_cv
load ff
load prc
load permno

fprintf('\n Now creating .csv for combination results. Run started at: ');
disp(datetime('now'));
tic;


% Get the anomalies; this step requires a lot of memory
[anoms,labels,anomaly_summary]=getChenZimmermanAnomalies();

% Pre-2005: 103; Continuous: 79; >50% in 1975: 58;
ind=anomaly_summary.Year<=2005 & strcmp(anomaly_summary.Cat_Form,'continuous');
chars=anoms(:,:,ind);
anomaly_summary2=anomaly_summary(ind,:);
labels2=labels(ind);

tempMe=me;
ind=makeUnivSortInd(me,5);
tempMe(ind==1)=nan;

indKeep=logical(zeros(size(labels')));

a=[];
nmonths=length(dates);

for i=1:size(chars,3)
    temp=chars(:,:,i);
    temp(isnan(tempMe))=nan;
    e=find(sum(isfinite(temp),2)>0,1,'last');
    if e>=nmonths-12 && e<nmonths
        temp(e:nmonths,:)=repmat(temp(e,:),nmonths-e+1,1);
    end
    temp2=sum(isfinite(temp),2)./sum(isfinite(tempMe),2);               
    s=find(temp2>=0.50,1,'first');
    if ~isempty(s) && s<=find(dates==197501)
        indKeep(i)=true;            
        temp(1:s,:)=nan;
    end     
    temp2=sum(isfinite(temp),2)./sum(isfinite(tempMe),2);                   
    a(:,i)=temp2;
    chars(:,:,i)=temp;                        
end

chars=chars(:,:,indKeep);
anomaly_summary2=anomaly_summary2(indKeep,:);
labels2=labels2(indKeep)';

% Convert to ranks
chars = tiedrank(permute(chars,[2 1 3]));
chars = (chars-1)./(max(chars)-1);
chars = chars-0.5;

% We are filling in the zeros & lagging here.
for i=1:size(chars,3)
    ind=find(sum(isfinite(chars(:,:,i)),1)>0);
    temp=chars(:,ind,i);
    temp(isnan(temp) & isfinite(tempMe(ind,:)'))=0;
    temp(:,1:find(dates(ind)==197501))=nan;
    chars(:,ind,i)=temp;
end

n=size(ret,1)*size(ret,2);
combination_data=[reshape(repmat(permno',size(ret,1),1),n,1) reshape(repmat(dates,1,size(ret,2)),n,1)];

for i=1:size(chars,3)
    combination_data=[combination_data reshape(chars(:,:,i)',n,1)];    
end

combination_data(sum(isfinite(combination_data),2)==2,:)=[];
combination_data=array2table(combination_data,'VariableNames',[{'permno','dates'}';labels2]);
writetable(combination_data,'Data/combination_data.csv');

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;

%%

clear
clc

load ret
load dates
load nyse
load me
load tcosts
load ff
load prc
load permno


fprintf('\n Now creating expected returns. Run started at: ');
disp(datetime('now'));
tic;

% Get the anomaly datasets
[anoms, labels]=getAnomalySignals('combination_data.csv','permno','dates');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the weighted-average signal expected returns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Now making the weighted-average signal expected returns at: ');
disp(datetime('now'));

% Run equal-weighted quintile sorts to determine the weights
for i=1:size(anoms,3)
    ind=makeUnivSortInd(anoms(:,:,i),5);
    res(i,1)=runUnivSort(ret,ind,dates,me,'weighting','e','timePeriod',197501,'printResults',0,'plotFigure',0);
end    
    
% Store the portfolio returns for the anomalies
s=find(dates==197501);
e=find(dates==202012);
pret=nan(length(dates),size(anoms,3));
for i=1:length(res)
    pret(s:e,i)=res(i).pret(:,end);    
end

eret=nan(size(ret));
for i=s+121:e
    weights=mean(pret(i-119:i,:),1)';
    weights(weights<0)=0;
    weights=weights/sum(weights);
    eret(i,:)=[permute(anoms(i,:,:),[2 3 1])*weights]';    
end

save Data/eret_waverage eret

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Fama-MacBeth expected returns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Now making the Fama-MacBeth expected returns at: ');
disp(datetime('now'));

eret=nan(size(ret));

T=120;

chars=permute(anoms,[2 1 3]);

for i=(T+1):length(dates)
    ind=find(sum(permute(sum(isfinite(chars(:,i-T+1:i,:)),1),[2 3 1])>0,1)==T);
    temp=permute(chars(:,i-T+1:i,ind),[2 1 3]);
    temp2=reshape(temp,size(temp,1),size(temp,2)*size(temp,3));
    res=runFamaMacBeth(100*ret(i-T+1:i,:),temp2,dates(i-T+1:i),'printResults',0);
    eret(i,:)=[permute(chars(:,i,ind),[1 3 2])*res.bhat(2:end)']';
end

save Data/eret_fmb eret

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Lasso expected returns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Now making the LASSO expected returns at: ');
disp(datetime('now'));

eret=nan(size(ret));

T=120;
chars=permute(anoms,[2 1 3]);

n=size(ret(find(dates==197501):end,:),1)*size(ret(find(dates==197501):end,:),2);

lassoY=reshape(ret(find(dates==197501):end,:),n,1);
lassoX=[];

for i=1:size(chars,3)
    temp=chars(:,find(dates==197412):end-1,i)';
    lassoX=[lassoX reshape(temp,n,1)];
end

lassoDates=reshape(repmat(dates(find(dates==197501):end),1,size(ret,2)),n,1);

ind=[isfinite(lassoY) & sum(isfinite(lassoX),2)>0];

lassoY=lassoY(ind);
lassoX=lassoX(ind,:);
lassoDates=lassoDates(ind);

bhat=nan(length(dates),size(chars,3));

for i=(find(dates==197501)+T):length(dates)
    tic;
    dates(i)
    ind=lassoDates<=dates(i) & lassoDates>=dates(i-119);
    y=lassoY(ind);
    x=lassoX(ind,:);
    [B,fitInfo]=lasso(x,y,'CV',5);
    minMSE=find(fitInfo.MSE==min(fitInfo.MSE));
    eret(i,:)=[permute(chars(:,i,:),[1 3 2])*B(:,minMSE)]';
    bhat(i,:)=B(:,minMSE)';
    toc;
end
save Data/eret_lasso eret bhat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the IPCA expected returns 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Now making the IPCA expected returns at: ');
disp(datetime('now'));

eret=nan(size(ret));

chars=anoms;

for i=1:size(chars,3)
    chars(:,:,i)=[lag(chars(:,:,i),1,nan)];
end

chars=permute(chars,[2 1 3]);

xret=(ret-repmat(rf,1,size(ret,2)))';

[~,~,L] = size(chars);
LOC = ~isnan(chars);
LOC = all(LOC,3)&~isnan(xret);
keepthese = sum(LOC)>=max(100,L+2);
LOC     = LOC(:,keepthese);
chars   = chars(:,keepthese,:);
xret    = xret(:,keepthese);
date    = dates(keepthese);

% Add a constant
chars(:,:,size(chars,3)+1)  = 1;
[N,T,L]         = size(chars);
Nts = sum(LOC);

% Construct X, W, X
chars   = permute(chars,[1 3 2]); % chars is now NxLxT
Z       = chars;
bigW       = nan(L,L,T);
bigX       = nan(L,T);
for t=1:T % parfor makes it too big
    bigW(:,:,t)    = (1/Nts(t)) * Z(LOC(:,t),:,t)' * Z(LOC(:,t),:,t); % W = Z'Z
    bigX(:,t)   = (1/Nts(t)) * Z(LOC(:,t),:,t)' * xret(LOC(:,t),t);  % X = Z'r
end


K=5;
als_opt.MaxIterations       = 10000;
als_opt.Tolerance           = 1e-6;
bigNts = sum(LOC);

[N,L,T]=size(Z);

temp_eret    = nan(N,T);

% Calculate restricted model first
% Initial guess for GammaBeta & the factor(s)    
for t=120:T-1
    X       = bigX(:,t-119:t);% this is X known through t
    W       = bigW(:,:,t-119:t);% this is W known through t
    Nts     = bigNts(t-119:t);% this is Nts known through t
    
    [GammaBeta_initial,s,v]    = svds(X,K);
    if t==120
        GB_Old      = GammaBeta_initial;
        F_Old       = s*v';%ones(K,T);%
    else
        GB_Old      = GB_New;
        F_Old       = [F_New(:,2:end) s*v(end,:)'];
    end
    
    % Iterate until you reach tolerance or max iterations
    tol         = 1;
    iter        = 0;
    tols        = nan(500,1);
    while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
        [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts);
        tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
        F_Old   = F_New;
        GB_Old  = GB_New;
        iter    = iter+1;
    end

    GammaBeta=GB_New;
    Lambda=mean(F_New,2);        
    
    temp_eret(:,t)     = Z(:,:,t+1)*GammaBeta*Lambda;    
end

eret(keepthese,:)=temp_eret';

save Data/eret_ipca eret

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;


%% Compare performance

clear
clc

load ret
load dates
load nyse
load me
load tcosts_cv
load ff

fprintf('\n Now comparing combination strategies. Run started at: ');
disp(datetime('now'));
tic;

% Make the low cost universe indicators
indME10=makeUnivSortInd(me,10,NYSE);
indLowCostTertile=nan(size(ret));
indLowCostHalf=nan(size(ret));

for i=1:10
    temp=tcosts;
    temp(indME10~=i)=nan;
    ind=makeUnivSortInd(temp,3);
    indLowCostTertile(ind==1)=1;
    ind=makeUnivSortInd(temp,2);
    indLowCostHalf(ind==1)=1;
end

JuneIndicator=(dates-100*floor(dates/100)==6);
JuneDecemberIndicator=(dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator=(dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

for j=1:4
    if j==1
        load eret_fmb
    elseif j==2
        load eret_ipca
    elseif j==3
        load eret_waverage
    elseif j==4
        load eret_lasso
    end
    for i=1:16
        if i==1 % Decile sort
            ind=makeUnivSortInd(eret,10);
        elseif i==2 % 10/20
            ind=makeUnivSortInd(eret,10,NYSE);
            ind=make_sS(ind,2);
        elseif i==3 % 10/30
            ind=makeUnivSortInd(eret,10,NYSE);
            ind=make_sS(ind,3);
        elseif i==4 % 10/40
            ind=makeUnivSortInd(eret,10,NYSE);
            ind=make_sS(ind,4);
        elseif i==5 %10/50
            ind=makeUnivSortInd(eret,10,NYSE);
            ind=make_sS(ind,5);
        elseif i==6 %20/25
            ind=makeUnivSortInd(eret,20);
            ind(ind==2 | ind==3 | ind==4)=1;
            ind(ind==19 | ind==18| ind==17)=20;
            ind=make_sS(ind,5);
        elseif i==7 % 20/30       
            ind = makeUnivSortInd(eret,10);    
            ind(ind==2)=1;
            ind(ind==9)=10;
            ind=make_sS(ind,3);
        elseif i==8 % 20/35
            ind = makeUnivSortInd(eret,20);    
            ind(ind==2 | ind==3 | ind==4)=1;
            ind(ind==19 | ind==18| ind==17)=20;
            ind=make_sS(ind,7);                
        elseif i==9 % 20/40
            ind = makeUnivSortInd(eret,10);    
            ind(ind==2)=1;
            ind(ind==9)=10;
            ind=make_sS(ind,4);               
        elseif i==10 % 20/45
            ind = makeUnivSortInd(eret,20);    
            ind(ind==2 | ind==3 | ind==4)=1;
            ind(ind==19 | ind==18| ind==17)=20;
            ind=make_sS(ind,9);                
        elseif i==11 % 20/50       
            ind = makeUnivSortInd(eret,10);    
            ind(ind==2)=1;
            ind(ind==9)=10;
            ind=make_sS(ind,5);      
        elseif i==12
            temp_eret=eret.*indLowCostHalf;
            ind=makeUnivSortInd(temp_eret,10);
        elseif i==13
            temp_eret=eret.*indLowCostTertile;
            ind=makeUnivSortInd(temp_eret,10);
        elseif i==14
            ind=makeUnivSortInd(eret,10);
            ind(~quarterEndIndicator,:)=0;
        elseif i==15
            ind=makeUnivSortInd(eret,10);
            ind(~JuneDecemberIndicator,:)=0;
        elseif i==16
            ind=makeUnivSortInd(eret,10);
            ind(~JuneIndicator,:)=0;
        end

        ind=1*(ind==1) + 2*(ind==max(max(ind)));
        resISew(j,i) = runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',[198501 200512],'plotFigure',0,'printResults',0); 
        resISvw(j,i) = runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',[198501 200512],'plotFigure',0);
        resOSew(j,i) = runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','e','timePeriod',[200512],'plotFigure',0,'printResults',0);
        resOSvw(j,i) = runUnivSort(ret,ind,dates,me,'tcosts',tcosts,'weighting','v','timePeriod',[200512],'plotFigure',0,'printResults',0);

    end
end

save Results/res_combo_strats resISew resISvw resOSew resOSvw

% Time-keeping

fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;

