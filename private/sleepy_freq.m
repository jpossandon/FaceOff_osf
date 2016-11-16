%   14 conditions
% 
% Correct
% 1 = 5_LT 6_RB w hair
% 2 = 5_RB 6_LT w hair
% 3 = 5_LB 6_RT w hair
% 4 = 5_RT 6_LB w hair
% 5 = 5_LT 6_RB w beard
% 6 = 5_RB 6_LT w beard
% 7 = 5_LB 6_RT w beard
% 8 = 5_RT 6_LB w beard
%

% 
%         efect/orthogonal coding
%         the model has 4 columns for position (center is coded -1)
%         1 column for face side (1/L, -1/R)
%         1 raw for freq (1/5hz, -1/6hz)
% 
%         
%         cont = [1 1 1 1 -1 -1 -1 -1 ;1 1 1 1 -1 -1 -1 -1];

%         pos = [1 -1 -1 1 1 -1 -1 1;-1 1 1 -1 -1 1 1 -1];

%         sid = [1 -1 1 -1 1 -1 1 -1;-1 1 -1 1 -1 1 -1 1];     

%% 

clear all;
close all; 

path        = '/Users/jossando/trabajo/FaceOff/';
experiment  = 'sleepy';
subj        = 1;
for reref       = 1;
    
    subj        = 1;
for session = {'S02_20131119','S03_20131120','S04_20131121','S05_20131122',...
         'S07_20131126','S08_20131126','S09_20131127','S10_20131128','S11_20131128',...
         'S12_20131129','S13_20131129','S15_20131203','S16_20131204','S17_20131204',...
         'S18_20131205','S19_20131209','S20_20131210',...
           'S21_20141014','S22_20141016','S23_20141016','S24_20141017','S25_20141021','S26_20141022','S27_20141023','S28_20141029','S29_20141030','S30_20141105'}%now the new batch,'S22_20141016}

    session     = session{:};
    EEGs        = pop_loadset('filename',[session '_Interpolated.set'],'filepath',[path 'data/' experiment '/' ]);
    stimulus    = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8'};
    
    
    % frequency analysis specs
      fs              = EEGs.srate;
        winsize         = 10; %(seconds)
        N               = winsize.*fs;
        cfg.method      = 'mtmfft';
        cfg.output      = 'pow';
        cfg.taper       = 'hanning';%'blackman'%
        cfg.keeptrials  = 'yes';
        cfg.foi         = 0:fs./N:fs.*(N-2)./(2.*N);
        cfg.foi         = cfg.foi(cfg.foi<25);
        cfg.t_ftimwin   = ones(1,length(cfg.foi))*winsize;
        cfg.toi         = 5;

     
    EEG             = pop_epoch( EEGs, stimulus(:), [-2  0], 'newname', 'ICA pruned with ICA epochs', 'epochinfo', 'yes');
    EEG             = pop_rmbase( EEG, [-1998 0]);
    
    EEG = pop_rmbase( EEG, [-1998 0]);
    if reref
        EEG = pop_reref( EEG, []);
    end
    databsl         = eeglab2fieldtrip( EEG, 'preprocessing', 'none');
    exlim           = floor(length(databsl.trial)*1000/5000);
    auxdata         = cell2mat(databsl.trial);
    databsl.trial   = mat2cell(auxdata(:,1:exlim*5000),64,5000*ones(1,exlim));
    databsl.time    = repmat({0:.002:9.998},1,exlim);
            [freq]                  = ft_freqanalysis(cfg, databsl);
    powbsl{subj}      = mean(freq.powspctrm);
    
    stimID      = 1;
    for stim = stimulus
        s   = stim{:};
        [~,rem] = strtok(s);
        numstim = str2num(rem);
       
        EEG = pop_epoch( EEGs, {  s  }, [0  10], 'newname', 'ICA pruned with ICA epochs', 'epochinfo', 'yes');

        EEG = pop_rmbase( EEG, [0  9998]);
        if reref
            EEG = pop_reref( EEG, []);
        end
         
        data            = eeglab2fieldtrip( EEG, 'preprocessing', 'none');

        [freq]                  = ft_freqanalysis(cfg, data);
        allpow(subj,stimID)       = freq;
        
        freqbl                  = freq;
        freqbl.powspctrm        = freq.powspctrm-repmat(powbsl{subj},[size(freq.powspctrm,1),1,1]);
        %baseline corrected absolute
        allpowbsl(subj,stimID)    = freqbl;
                 
        freq.powspctrm          = freq.powspctrm./repmat(powbsl{subj},[size(freq.powspctrm,1),1,1]);
        %baseline corrected relative change
        allpowbslrl(subj,stimID)    = freq;
        
        stimID          = stimID+1;
    end
   subj = subj+1; 
 
end

    if reref
        save([path 'results/allfreqsnew_080816_' experiment '_AR'],'allpow','allpowbsl','allpowbslrl')
    else
        save([path 'results/allfreqsnew_080816_' experiment '_NR'],'allpow','allpowbsl','allpowbslrl')
    end
      clear allpow allpowbsl allpowbslrl
end


  
    
  