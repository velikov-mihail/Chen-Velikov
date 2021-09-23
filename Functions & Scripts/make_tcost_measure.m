clear
clc

fprintf('\n Now creating trading cost measures. Run started at: ');
disp(datetime('now'));
tic;


gibbsFileName='crspgibbs2021.csv';

effSpreadStruct=struct;
effSpreadStruct.gibbs=makeGibbs(gibbsFileName); % Returns 2*c (the full spread). We'll divide by 2 later. Should have crspgibbsYYYY.csv somewhere on the path.
effSpreadStruct.hl=makeCorwinSchultz();
effSpreadStruct.chl=makeAbdiRanaldi();
effSpreadStruct.vov=makeKyleObizhaeva();
effSpreadStruct.hf_spreads=makeHighFreqEffSpreads();
vars={'gibbs','chl','hl','vov','hf_spreads.ave','hf_spreads.monthend'};

save Data/effSpreadStruct_cv effSpreadStruct -v7.3 % Adding _cv at the end to distinguish from tcost variables stored by the MATLAB Asset Pricing library

% Store some dimensions
nummonths=size(effSpreadStruct.gibbs,1);
numstocks=size(effSpreadStruct.gibbs,2);
numobs=nummonths*numstocks;

% Winsorize all measures at the top at 99.9%
for i=1:length(vars)  
    eval(['temp=effSpreadStruct.',char(vars(i)),';']);
    temp2=reshape(temp,numobs,1);
    ind=isfinite(temp2);
    temp2(ind)=winsor(temp2(ind),[0 99.9]);    
    temp3=reshape(temp2,nummonths,numstocks);
    eval(['effSpreadStruct.',char(vars(i)),'=temp3;']);    
end

load dates
load exchcd

% Exlude NASDAQ stocks prior to 1983 for Gibbs and VoV
ind=(repmat(dates,1,size(exchcd,2))<=1982212 & exchcd==3);
effSpreadStruct.gibbs(ind)=nan;
effSpreadStruct.vov(ind)=nan;

% Exlude NASDAQ stocks prior to 1993 for HL and CHL
ind=(repmat(dates,1,size(exchcd,2))<=199212 & exchcd==3);
effSpreadStruct.hl(ind)=nan;
effSpreadStruct.chl(ind)=nan;

% Exlude AMEX stocks prior to 1962 for all
ind=(repmat(dates,1,size(exchcd,2))<=196112 & exchcd==2);
effSpreadStruct.gibbs(ind)=nan;
effSpreadStruct.hl(ind)=nan;
effSpreadStruct.chl(ind)=nan;
effSpreadStruct.vov(ind)=nan;        

% Combine the measures
reshapedLF=[reshape(effSpreadStruct.gibbs,numobs,1) reshape(effSpreadStruct.hl,numobs,1) ...
            reshape(effSpreadStruct.chl,numobs,1) reshape(effSpreadStruct.vov,numobs,1)];
reshapedEffSpreadRaw=nanmean(reshapedLF,2); % Use nanmean here, so that if at least one non-NaN we'll use it
reshapedHF=reshape(effSpreadStruct.hf_spreads.ave,numobs,1);        
reshapedEffSpreadRaw(isfinite(reshapedHF))=reshapedHF(isfinite(reshapedHF)); % Assign high-frequency spreads wherever we have it        
effSpreadRaw=reshape(reshapedEffSpreadRaw,nummonths,numstocks); % Turn back into a matrix 

tcosts_raw=effSpreadRaw/2; % Need to divide by 2, because this is the tcost measure (half-spread!!!)
save Data/tcosts_raw_cv tcosts_raw

tcosts=fillMissingTcosts(tcosts_raw);
save Data/tcosts_cv tcosts


% Time-keeping
fprintf('\n Run ended at: ');
disp(datetime('now'));
toc;
