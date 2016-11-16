%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for the linear model analysis of
% osf.io project: Face Off
% This code reproduces the results presented on the article:
%
% Requires the re-processed spectral data obtained with script:
%   sleepy_freq.m 
% from each subject eeglab sets (which have been visually and ICA cleaned).
% The resulting file is named:
%       allfreqsnew_<date>_sleepy_AR.mat
% and contain three structure arrays (allpow,allpowbsl,allpowbslrl). Every
% structure array has (#subjects,#condition) components which are the result 
% of spectral analysis performed with functions from the fieldtrip toolbox 
% (as specified in sleepy_freq.m script). 
% 
% The 8 conditions are:
%       condition = 5Hz face side (Left/Right) and position (Top/Bottom) 
%                   / same 6 Hz / facial contrast (hair or beard)
%       
%       1 = 5_LT 6_RB w hair
%       2 = 5_RB 6_LT w hair
%       3 = 5_LB 6_RT w hair
%       4 = 5_RT 6_LB w hair
%       5 = 5_LT 6_RB w beard
%       6 = 5_RB 6_LT w beard
%       7 = 5_LB 6_RT w beard
%       8 = 5_RT 6_LB w beard
% 
% The mapping from the structure arrays to the model is descrived in three 
% varialbles (cont,pos,side) whithin the strcuture coding. Columns correspond
% to the conditions described above and the rows to the two different 
% frequencies:
% 
%     .cont = [1 1 1 1 -1 -1 -1 -1 ;...          % facial contrast
%              1 1 1 1 -1 -1 -1 -1];             % hair (1) or beard (-1)
% 
%     .pos = [1 -1 -1  1  1 -1 -1  1;...         % position
%            -1  1  1 -1 -1  1  1 -1];           % top (1) or bottom (-1)
% 
%     .sid = [1 -1  1 -1  1 -1  1 -1;...         % face side
%            -1  1 -1  1 -1  1 -1  1];           % left (1) or right (-1)
%
%     The glm analysis if done with efect/orthogonal coding, 
%      the model design matrix has:
% 
%         1 column for position         (pos : 1/top  , -1/bottom)
%         1 column for face side        (sid : 1/L    , -1/R)
%         1 column for frequency        (freq: 1/5hz  , -1/6hz)
%         1 column for facial contrast  (cont: 1/hair , -1/beard)
%         6 columns for the two-way interactions
%         the model also generates a result for a constant (grand mean), but
%         this is not entered as a column in the desgin matrix variable (XY)
% 
% by Jose Ossandon (jose.ossandon@uni-hamburg.de)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
% Linear model

clear 
% coding
coding.cont     = [1 1 1 1 -1 -1 -1 -1 ;1 1 1 1 -1 -1 -1 -1];
coding.pos      =  [1 -1 -1 1 1 -1 -1 1;-1 1 1 -1 -1 1 1 -1];
coding.sid      = [1 -1 1 -1 1 -1 1 -1;-1 1 -1 1 -1 1 -1 1];         

% used frequencies
 fr             = [5;6];  % 5/6 Hz  
         
path            = '/Users/jossando/trabajo/FaceOff/code/osf/';
experiment      = 'sleepy';

load([path 'results/allfreqsnew_080816_' experiment '_AR'])

times           = allpow(1,1).cfg.toi;
frequencies     = allpow(1,1).freq;
datestr         = datestr(now,'ddmmyy');
npermute        = 0;                                                        % within subject permutations
npermuteB       = 10000;                                                    % between subjects permutation
load([path 'channels/channel_loc']) 
     
%%
statname = 'signpermT';                                                     %''bootet';%bootet';%'bootet';%'signpermT';
mcname   = 'cluster';                                                       %'tfce';%;'maxsT';%
dataname = 'allpow';                                                        % spectral analysis and then average
auxdata  = allpow;
   
Y = {}; XY = {};
for suj = 1:size(auxdata,1)
    side = []; cont = []; pos = []; frequ = []; Yaux = [];
    for f = 1:size(fr,1)
        for numstim = 1:size(auxdata,2)
            % get data for a given subject, trigger, and frequency
            indxfreq            = find(frequencies==fr(f));
            dump                = squeeze(sum(auxdata(suj,numstim).powspctrm(:,:,indxfreq,:),3));

            nT(suj,numstim,f)   = size(dump,1);    % amount of trials per condition and subject
            if isempty(Yaux)       % dependent variable              
                Yaux            = dump;
            else
                Yaux            = cat(1,Yaux,dump); 
            end
            % effect coding
            side                = [side;coding.sid(f,numstim)*ones(nT(suj,numstim,f),1)];
            cont                = [cont;coding.cont(f,numstim)*ones(nT(suj,numstim,f),1)];
            pos                 = [pos;coding.pos(f,numstim)*ones(nT(suj,numstim,f),1)];
            if f == 1
                frequ          = [frequ;ones(nT(suj,numstim,f),1)];
            elseif f == 2
                frequ          = [frequ;-ones(nT(suj,numstim,f),1)];     
            end

        end
    end
    XY{suj} = [side,cont,pos,frequ,side.*cont,side.*pos,side.*frequ,cont.*pos,cont.*frequ,pos.*frequ];
    Y{suj}  = Yaux;
end

%%
for np = 1:npermute+1
    tstart = tic;
    for suj = 1:size(auxdata,1)

        % subject linear model
        XYaux = XY{suj};
        if np>1
            XYaux = XYaux(randsample(1:size(XYaux,1),size(XYaux,1)),:);
        end
        [modelaux(suj).B,modelaux(suj).Bt,modelaux(suj).STATS,model(suj).TCFE] = regntcfe(Y{suj},XYaux,1,'effect',elec,1);

        % side
        [~,~,models.sid(suj).STATS]     = regntcfe(Y{suj},XYaux(:,1),1,'effect',elec,1);
        % facial contrast
        [~,~,models.cont(suj).STATS]    = regntcfe(Y{suj},XYaux(:,2),1,'effect',elec,1);
        % position
        [~,~,models.pos(suj).STATS]     = regntcfe(Y{suj},XYaux(:,3),1,'effect',elec,1);
        % freq
        [~,~,models.freq(suj).STATS]    = regntcfe(Y{suj},XYaux(:,4),1,'effect',elec,1);
        % interactions
        [~,~,models.inter(suj).STATS]   = regntcfe(Y{suj},XYaux(:,5:end),1,'effect',elec,1);
        
        if suj ==1
            betasaux = modelaux(suj).B;
        else
            betasaux = cat(4,betasaux,modelaux(suj).B);
        end

    end
    
    if np==1
        [result]    = regmodel2ndstat(betasaux,times,elec,npermuteB,statname,mcname);
        model       = modelaux;
        betas       = betasaux;
        if ~strcmp(statname,'signpermT')
            for b = 1:size(betasaux,2)
            [clustersac]        = clustereeg(reshape(result.T(:,:,b),[1,64]),reshape(result.Hnc(:,b),[1,64]),elec,64,1);
            result.clusters(b)  = clustersac;
            end
        end
    else
        [npresult]  = regmodel2ndstat(betasaux,times,elec,0,statname,mcname);

        for b = 1:size(betasaux,2)
        [clustersnp]            = clustereeg(reshape(npresult.T(:,:,b),[1,64]),reshape(npresult.Hnc(:,b),[1,64]),elec,64,1);
         npMAXclus(np-1,b)      = clustersnp.MAXst;
        end
    end
    sprintf('%s Permutation %d model %s %4.2f',experiment,np,dataname,toc(tstart))
end

%%
if np>1
    for b = 1:size(betas,2)
         result.clusters(b).npMAXclus =  npMAXclus(:,b);
    end
end
save([path 'results/glm_AR_' experiment '_' datestr '_' dataname '_' statname '_' mcname ],'result','betas','model','coding','nT','models')

