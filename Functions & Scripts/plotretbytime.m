function [lez,lsmooth,lse1] = plotretbytime(time,ret,Tsmooth)
% 2017 11

% param
SEfac = 2;
% ret(:,[33 38 141])=[];
% calculate stats
retez = meanez(ret,2);
[Eretsmooth,SEsmooth,Nport] = nanall(size(ret,1),1);
for t = Tsmooth+1:length(ret)
    tempret = ret(t-Tsmooth:t,:);
    tempvec = cropnan(tempret(:));
    Eretsmooth(t) = mean(tempvec(:));        
    SEsmooth(t) = std(tempvec)/sqrt(length(tempvec));
%         Nport(t) = length(tempvec)/Tsmooth;
end    

retez=nanmean(ret,2);
Eretsmooth=moving(retez,60);


% plot
lez = plot(time,retez);       
temp = time(1):time(end);
l0 = plot(temp,zeros(size(temp)),'k-'); % plot zero line    
vline(0,'k-')
lsmooth = plot(time,Eretsmooth);
lse1 = plot(time,Eretsmooth+SEfac*SEsmooth,'r--');
lse2 = plot(time,Eretsmooth-SEfac*SEsmooth,'r--');    

set(lez     ,'color',[1 1 1]*0.8)
set(lsmooth ,'color',[   0    0.4470    0.7410],...
    'linewidth',1.5)
set(lse1     ,'color',[0.9 0 0])
set(lse2     ,'color',[0.9 0 0])    


