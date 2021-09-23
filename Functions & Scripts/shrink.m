function [point,boot] = shrink(samp,numu,Nboot,plist)


% - SR is t-dist
sr = samp.SR';
sesr = 1./sqrt((samp.tstat'./samp.SR').^2);
vol = samp.vol';


ok = ~isnan(samp.SR);
sr = sr(ok); sesr = sesr(ok); vol = vol(ok);

N = length(sr);


rng(1)
anomid = randi([1 N], [N Nboot]);

[boot.radj_quant, boot.radj_quant_vol] = nanall(Nboot,4);
for booti = 0:Nboot
    
    % draw anomalies
    if booti == 0
        % here don't randomize
        bootid = 1:N;
    else
        bootid = anomid(:,booti);
    end
    
    bootsr = sr(bootid);
    bootsesr = sesr(bootid);
    bootvol = vol(bootid);
    
    % estimate
    musr = mean(bootsr);
    tempv = (numu-2)/numu*(var(bootsr)-mean(bootsesr.^2));
    tempv = max(tempv,0);
    sigsr = sqrt(tempv);    
    
    % shrink
    [bootsr_adj, muhatcheck] = nanall(size(bootsr));
    if sigsr > 0
        for i = 1:length(bootsr)
            
            curr_mean = bootsr(i); curr_sig = bootsesr(i);
            
            kern_denom = @(mu) normpdf(curr_mean,mu,curr_sig) ...
                .*(1/sigsr).*tpdf((mu-musr)./sigsr, numu);
            
            kern_num = @(mu) mu.*kern_denom(mu);
            
            bootsr_adj(i) = integral(kern_num,-inf,inf)/integral(kern_denom,-inf,inf);
            
            % check
            s = curr_sig^2/(curr_sig^2+sigsr^2);
            muhatcheck(i) = s*musr+(1-s)*curr_mean;
            
        end
    else
        % 100 % shrinkage
        bootsr_adj = musr;
    end
    
    % find adjusted returns
    bootr_adj = bootsr_adj.*bootvol;
    
    % summarize
    bootradj_quant = ptile(bootr_adj, plist);   

    % look up vol corresponding to these bootr_adj
    dist = abs(bsxfun(@plus, bootr_adj',  -bootradj_quant));
    [~,istar] = min(dist);
    bootradj_quant_vol = bootvol(istar);      

    
    % store
    if booti == 0
        point.radj_quant = bootradj_quant;
        point.radj_quant_vol = bootradj_quant_vol;
        point.sigsr = sigsr;
        point.musr = musr;
        
    else
        boot.radj_quant(booti,:) = bootradj_quant;
        boot.radj_quant_vol(booti,:) = bootradj_quant_vol;
        boot.sigsr(booti,1) = sigsr;
        boot.musr(booti,1) = musr;        
    end
    
    
end % for booti
