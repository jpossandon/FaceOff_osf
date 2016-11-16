%   14 conditions
% 
%          5Hz left / 6Hz right
%         1 = -3M -4L -2R
%         2 = -2M -3L -1R
%         3 = -1M -2L  0R
%         4 =  0M -1L  1R
%         5 =  1M  0L  2R
%         6 =  2M  1L  3R
%         7 =  3M  2L   4R
%  
%          6Hz left / 5Hz right
%         8 =  -3M -4L -2R
%         9 =  -2M -3L -1R
%         10 = -1M -2L  0R
%         11 =  0M -1L  1R
%         12 =  1M  0L  2R
%         13 =  2M  1L  3R
%         14 =  3M  2L  4R
% 
%         efect coding
%         the model has 4 columns for position (center is coded -1)
%         1 column for face side (1/L, -1/R)
%         1 colun for freq (1/5hz, -1/6hz)
% 
%         
%         dat = [0 0 1 1 1 1 1 1 1 1 1 1 0 0 ; -1 -1 -1 -1 -1 0 0 0 0 -1 -1 -1 -1 -1];
%         pos = [0 0 1 2 5 3 4 1 2 5 3 4 0 0;1 2 5 3 4 0 0 0 0 1 2 5 3 4];
%         sid = [1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1;-1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1];     

%% 
% TODO: baseline
 clear;
close all; 

path        = '/Users/jossando/trabajo/FaceOff/';
experiment  = 'eccentricity';

for reref       = 1;
subj        = 1;
for session = {'S02_20131119','S03_20131120','S04_20131121','S05_20131122',...
         'S07_20131126','S08_20131126','S09_20131127','S10_20131128','S11_20131128',...
         'S12_20131129','S13_20131129','S15_20131203','S16_20131204','S17_20131204',...
         'S18_20131205','S19_20131209','S20_20131210',... 
         'S21_20141014','S22_20141016','S23_20141016','S24_20141017','S25_20141021','S26_20141022','S27_20141023','S28_20141029','S29_20141030','S30_20141105'}  %these are eeglab set of continuous data visually and ICA cleaned

    session     = session{:};
    EEGs = pop_loadset('filename',[session '_Interpolated.set'],'filepath',[path 'data/' experiment '/' ]);
        
    stimulus    = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8','S  9','S 10','S 11','S 12','S 13','S 14';...
                    'S 1','S  2','S 3','S  4','S 5','S  6', 'S 7','S  8','S 9','S 10','S 11','S 12','S 13','S 14'};
    stimID      = 1;
    
    % freq analysis specs
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
        
    %baseline calculation
%=     if strcmp(session,'S22_20141016')
%         EEG             = pop_epoch( EEGs, stimulus(2,:), [-2  0], 'newname', 'ICA pruned with ICA epochs', 'epochinfo', 'yes');
%     else
         EEG             = pop_epoch( EEGs, stimulus(1,:), [-2  0], 'newname', 'ICA pruned with ICA epochs', 'epochinfo', 'yes');
%     end
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
    
    for stim = 1:size(stimulus,2)
%         if strcmp(session,'S22_20141016')
%             s   = stimulus{2,stim};
%         else
            s   = stimulus{1,stim};
%         end
        [~,rem] = strtok(s);
        numstim = str2num(rem);
        

            
        EEG = pop_epoch( EEGs, {  s  }, [0  10], 'newname', 'ICA pruned with ICA epochs', 'epochinfo', 'yes');
%         end
        EEG = pop_rmbase( EEG, [0  9998]);
        if reref
            EEG = pop_reref( EEG, []);
        end
        data            = eeglab2fieldtrip( EEG, 'preprocessing', 'none');

      

        [freq]                  = ft_freqanalysis(cfg, data);
        allpow(subj,stim)       = freq;
        
        freqbl                  = freq;
        freqbl.powspctrm          = freq.powspctrm-repmat(powbsl{subj},[size(freq.powspctrm,1),1,1]);
        %baseline corrected
        allpowbsl(subj,stim)    = freqbl;
                 
        freq.powspctrm          = freq.powspctrm./repmat(powbsl{subj},[size(freq.powspctrm,1),1,1]);
        %baseline corrected  10*log10(data ./ meanVals)
        allpowbslrl(subj,stim)    = freq;
       
        
         datatl                  = ft_timelockanalysis([],data);
         [freqtl]                = ft_freqanalysis(cfg, datatl);
         allpowtl(subj,stim)   = freqtl;
       
%         stimID          = stimID+1;
    end
    subj = subj+1;
    
end

    if reref
        save([path 'results/allfreqsnew_260416_' experiment '_AR'],'allpow','allpowbsl','allpowbslrl','allpowtl')
    else
        save([path 'results/allfreqsnew_260416_' experiment '_NR'],'allpow','allpowbsl','allpowbslrl','allpowtl')
    end
    clear allpow allpowbsl allpowbslrl
end

  
    
  