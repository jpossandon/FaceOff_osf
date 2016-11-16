function [B,Bt,STATS,T] = regntcfe(Y,XY,p,coding,elec,subjst)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B,Bt,STATS,T] = regntcfe(Y,XY)
% Does reression of every column of Y according to model XY (without
% constant, and then 2d tcfe of the t values of each coficients according
% to the neighbour structure within elec
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B       = NaN([size(Y,2),size(XY,2)+1,size(Y,3)]);  
Bt      = NaN([size(Y,2),size(XY,2)+1,size(Y,3)]);
STATS   = NaN([size(Y,2),7,size(Y,3)]);
mask    = randsample([-1 1],size(Y,1),'true')';  % this was before inside the loop which was wrong because it decorrelated the temporal structure

for ch = 1:size(Y,2)
%     sprintf('Processing Channel %d/%d',ch,size(Y,2))
    for t = 1:size(Y,3)
        stats           = regstats(Y(:,ch,t),XY,'linear',{'tstat','fstat','rsquare'});
        Bt(ch,:,t)      = stats.tstat.t;
        B(ch,:,t)      = stats.tstat.beta;
        STATS(ch,:,t)  = [stats.rsquare,stats.fstat.f,stats.fstat.pval,stats.fstat.sse,stats.fstat.dfe,stats.fstat.ssr,stats.fstat.dfr];
        

        if p>1 && strcmp(coding,'effect')% to test the constant with effect coding
            stats       = regstats(Y(:,ch,t).*mask,XY,'linear',{'tstat'});
            Bt(ch,1,t)  = stats.tstat.t(1);
        end
    end
    
end
if subjst
    for b = 1:size(Bt,2) % This should be used only for single subject statistics
        if size(Y,2)>1
            T(:,:,b) = tfce(squeeze(Bt(:,b,:)),[],elec.channeighbstructmat,'stat');
        else
            T(:,:,b) = tfce(squeeze(Bt(:,b,:))',[],elec.channeighbstructmat,'stat');
        end
    end
else
    T=[];
end
