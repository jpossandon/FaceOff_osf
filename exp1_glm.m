%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for the linear model analysis of
% osf.io project: Face Off
% This code reproduces the results presented on the article:
%
% Requires the re-processed spectral data obtained with script: 
%   eccentricity_freq.m 
% from each subject eeglab sets (which have been visually and ICA cleaned).
% The resulting file is named:
%       allfreqsnew_<date>_eccentricity_AR.mat
% and contain three structure arrays (allpow,allpowbsl,allpowbslrl). Every
% structure array has (#subjects,#condition) components which are the result 
% of spectral analysis performed with functions from the fieldtrip toolbox 
% (as specified in eccentricity_freq.m script). 
% 
% The 14 conditions are:
%      condition = whole head position / left hemiface position
%                 / right hemiface pos)
%         
%           5Hz left / 6Hz right
%                 1 = -3M -4L -2R
%                 2 = -2M -3L -1R
%                 3 = -1M -2L  0R
%                 4 =  0M -1L  1R
%                 5 =  1M  0L  2R
%                 6 =  2M  1L  3R
%                 7 =  3M  2L   4R
%  
%           6Hz left / 5Hz right
%                 8  = -3M -4L -2R
%                 9  = -2M -3L -1R
%                 10 = -1M -2L  0R
%                 11 =  0M -1L  1R
%                 12 =  1M  0L  2R
%                 13 =  2M  1L  3R
%                 14 =  3M  2L  4R
% 
% The mapping from the structure arrays to the model is descrived in three 
% varialbles (dat,pos,side), in which columns correspond to the conditions 
% described above and the rows to the two different frequencies:
%
%   .dat = [0  0  1  1  1 1 1 1 1  1  1  1  0  0;...     % whereas there is data for the given position (column) 
%          -1 -1 -1 -1 -1 0 0 0 0 -1 -1 -1 -1 -1];       % and frequency (row1:5Hz, row2:6Hz)
%
%   .pos = [0 0 1 2 5 3 4 1 2 5 3 4 0 0;...              % which position, with respect to the coding described above 
%           1 2 5 3 4 0 0 0 0 1 2 5 3 4];                % the mapping is (-2:1, -1:2, 0:5, 1:3, 2:4)
%
%   .sid = [1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1;... % face side
%          -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1];   % left (1) or right (-1)
%
%     The glm analysis if done with efect/orthogonal coding, 
%      the model design matrix has:
%         4 columns for position (1st column:-2, 2nd:-1, 3rd:1, 4th:2, and center is coded with -1 in all columns)
%         1 column for face side (1/L, -1/R)
%         1 colun for freq (1/5hz, -1/6hz)
%         9 columns for the two-way interactions
%         the model also generates a result for a constant (grand mean), but
%         this is not entered as a column in the design matrix variable (XY)
% 
% by Jose Ossandon (jose.ossandon@uni-hamburg.de) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%
% Linear model

clear 
% coding
coding.dat = [0 0 1 1 1 1 1 1 1 1 1 1 0 0 ; -1 -1 -1 -1 -1 0 0 0 0 -1 -1 -1 -1 -1];
coding.pos = [0 0 1 2 5 3 4 1 2 5 3 4 0 0;1 2 5 3 4 0 0 0 0 1 2 5 3 4];
coding.sid = [1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1];     


fr              = [5;6];                                                    % used frequencies
path            = '/Users/jossando/trabajo/FaceOff/code/osf/';
experiment      = 'eccentricity';
datestr         = datestr(now,'ddmmyy');
npermute        = 0;                                                        % single subject permutations, set to 1 for main analysis
npermuteB       = 10000;                                                    % between subjects permutation

load([path 'results/allfreqsnew_260416_eccentricity_AR'])                   % this is the pre-processed spectral data obtained with script eccentricity_freq.m from eeglab set (which have been visually and ICA cleaned)
times           = allpow(1,1).cfg.toi;
frequencies     = allpow(1,1).freq;
load([path 'channels/channel_loc']) 
     
%%
% getting different subjects data
statname    = 'signpermT';                                                  % 2nd level analysis is done through permutation of the sign of subject's betas
mcname      = 'cluster';                                                    % values are clustered across space
dataname    = 'allpow';                                                     % spectral analysis and then average
auxdata     = allpow;
 
Y = {}; XY = {};
for suj = 1:size(auxdata,1)
    side = []; posc = []; frequ = []; Yaux = [];
    for f = 1:size(fr,1)
        for numstim = 1:size(auxdata,2)

            % get data for a given subject, trigger, and frequency
            indxfreq            = find(ismember(frequencies,fr(f,:)));
            dump                = squeeze(sum(auxdata(suj,numstim).powspctrm(:,:,indxfreq,:),3));
            
            if coding.dat(f,numstim)
                nT(suj,numstim,f)   = size(dump,1);    % amount of trials per condition and subject
                if isempty(Yaux)       % dependent variable              
                    Yaux       = dump;
                else
                    Yaux       = cat(1,Yaux,dump); 
                end
                % effect coding
                side        = [side;coding.sid(f,numstim)*ones(nT(suj,numstim,f),1)];
                posaux      = zeros(nT(suj,numstim,f),4);
                if coding.pos(f,numstim)>0 && coding.pos(f,numstim)<5
                    posaux(:,coding.pos(f,numstim)) = 1;
                elseif coding.pos(f,numstim)==5
                    posaux(:)                       = -1;
                end
                posc        = [posc;posaux];
                if f == 1
                    frequ   = [frequ;ones(nT(suj,numstim,f),1)];
                elseif f == 2
                    frequ   = [frequ;-ones(nT(suj,numstim,f),1)];     
                end
            end
        end
    end
    
    XY{suj}     = [side,posc,frequ,repmat(side,[1,4]).*posc,...
                    side.*frequ,repmat(frequ,[1,4]).*posc];                 % predictors
    Y{suj}      = Yaux;                                                     % data
        
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
        [modelaux(suj).B,modelaux(suj).Bt,modelaux(suj).STATS,modelaux(suj).TCFE] = regntcfe(Y{suj},XYaux,1,'effect',elec,1);
        
        % reduced models (orthogonal design, we can evalute every factor independently)
        % side
        [~,~,models.sid(suj).STATS] = regntcfe(Y{suj},XYaux(:,1),1,'effect',elec,1);
        %pos
        [~,~,models.pos(suj).STATS] = regntcfe(Y{suj},XYaux(:,2:5),1,'effect',elec,1);
        %freq
        [~,~,models.freq(suj).STATS] = regntcfe(Y{suj},XYaux(:,6),1,'effect',elec,1);
        %interact
        [~,~,models.inter(suj).STATS] = regntcfe(Y{suj},XYaux(:,7:end),1,'effect',elec,1);
        %check sidexpos
        [~,~,models.sidepos(suj).STATS] = regntcfe(Y{suj},XYaux(:,1:5),1,'effect',elec,1);
        
        % calculate the position in the middle
        if suj ==1
            betasaux = [modelaux(suj).B(:,1:4),-1*sum(modelaux(suj).B(:,3:6),2),...
                modelaux(suj).B(:,5:9),-1*sum(modelaux(suj).B(:,8:11),2),...
                modelaux(suj).B(:,10:14),-1*sum(modelaux(suj).B(:,13:16),2),...
                modelaux(suj).B(:,15:16)];
        else
            betasaux = cat(4,betasaux, [modelaux(suj).B(:,1:4),...
                -1*sum(modelaux(suj).B(:,3:6),2),modelaux(suj).B(:,5:9),...
                -1*sum(modelaux(suj).B(:,8:11),2),modelaux(suj).B(:,10:14),...
                -1*sum(modelaux(suj).B(:,13:16),2),modelaux(suj).B(:,15:16)]);
        end
             
    end
    % 2nd level analysis
    if np==1
        model       = modelaux;
        betas       = betasaux;
        [result]    = regmodel2ndstat(betasaux,times,elec,npermuteB,statname,mcname);
        if ~strcmp(statname,'signpermT')
            for b = 1:size(betasaux,2)
            [clustersac]        = clustereeg(reshape(result.T(:,:,b),[1,64]),reshape(result.Hnc(:,b),[1,64]),elec,64,1);
            result.clusters(b)  = clustersac;
            end
        end
    else
        [npresult]  = regmodel2ndstat(betasaux,times,elec,0,statname,mcname);
        for b = 1:size(betasaux,2)
           [clustersnp]         = clustereeg(reshape(npresult.T(:,:,b),[1,64]),reshape(npresult.Hnc(:,b),[1,64]),elec,64,1);
            npMAXclus(np-1,b)   = clustersnp.MAXst;
        end
    end
        sprintf('%s Permutation %d model %s %4.2f',experiment,np,dataname,toc(tstart))
end
if np>1
    for b = 1:size(betas,2)
         result.clusters(b).npMAXclus =  npMAXclus(:,b);
    end
end
save([path 'results/glm_AR_' experiment '_' datestr '_' dataname '_' statname '_' mcname ],'result','betas','model','coding','nT','models')

