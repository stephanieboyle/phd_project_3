%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING, behaviour, ERPS, DECODING, stimulus locked analysis
% PREPROCESSING, DECODING, response locked analysis
% steph 3/11/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
cd('W:\CIRCLE_EXP\LOG\EEG');
addpath('W:\CIRCLE_EXP');
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')

% subs = dir('Subj_*');

%--------------------------------------------------------------------------
% PREPRO PARAMETERS
ARG.Amp_thr = 120;                      % threshold for complete trial rejection
ARG.Resample = 200;                     % sampling rate, use 150 or 200
ARG.Highpass = 0.5;                     % high pass filter

% Artifact Removal parameters
ARG.Remove = 'all';                     % 'fast' for  fast processing and removal of eye topo only 'all' for removal of other noise sources
ARG.ICA_artif_thr = 0.7;                % threshold: correlation with template
ARG.ICA_eyecorr = 0.08;                 % threshold : correlation with eye movement signal (only for 'all')
ARG.ICA_spectrum_thr = 6;               % spectrum on component amplitude (only for 'all')


for SUBJ = input('Subj: \n') %1:length(subs);
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d/INDIVIDUAL_BLOCKS',SUBJ))
    Blocks = dir('CIRCLE_S*_B*_*.mat');
    
    for B = 1:length(Blocks);  % BLOCK
        
        %%
        
        B = input('Block? \n');
        fprintf('\n\n\n')
        fprintf('Subj %d, Block %d... \n',SUBJ,B);
        fnameB = Blocks(B).name;
        fnameEEG = dir(sprintf('CIRCLE_S%d_B%d_*.bdf',SUBJ,B));
        
        % TRIAL PARAMETERS
        trialdef.prestim = -0.5;                                                % pre trigger secs
        trialdef.poststim = 1.5;                                                % after trigger secs
        trigRange = 1;
        trigs = length(trigRange);                                              % for later (to work out if there are 2 sets of triggers or 1)
        
        % REREFERENCING
        reref = 0;                                                              % none
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % READ EVENTS AND HEADER FILES
        fullname = fnameEEG.name;                                               % EEG FILE
        matname = fnameB;                                                       % BEHAVIOURAL MATLAB FILE
        cfg = [];
        cfg.dataset = fullname;                                                 % EEG data.
        event = ft_read_event(cfg.dataset);
        hdr   = ft_read_header(cfg.dataset);
        EVENTS = ft_filter_event(event, 'type','STATUS');                       % find all triggers
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % FIND TRIGGERS
        clear val t_onset
        c=1;
        %         val = zeros(1,160); t_onset = val;                                      % initialise
        for k=1:length(EVENTS)
            if sum(EVENTS(k).value==trigRange) && (EVENTS(k).sample>1)
                val(c) = EVENTS(k).value;                                       % trigger value
                t_onset(c) = EVENTS(k).sample;                                  % sample number
                c=c+1;
            end
        end
        fprintf('Found %d trials \n',c-1);
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % LOAD IN MATLAB DATA
        load(matname)
        A = load(matname,'Response');
        %         A = load(matname,'all');
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % PUT DATA INTO FIELDTRIP TRIAL STRUCTURE
        %         trl =zeros(size(Response,1),size(Response,2)+4);
        trl = [] ;
        for t=1:length(t_onset)
            ts =t_onset(t);
            begsample     = ts + trialdef.prestim*hdr.Fs;                       % sample at trigger
            endsample     = ts + trialdef.poststim*hdr.Fs;                      % sample at end
            offset        = trialdef.prestim*hdr.Fs;                            % offset
            trl(t,(1:3)) = round([begsample endsample offset]);
            trl(t,4) = val(t);                                                  % add in the trigger values
        end
        trl = cat(2,trl,Response);
        %         trl = cat(2,trl,all);
        %------------------------------------------------------------------
        
        % subject 10 was so fast there wasn't enough data on the last trial
        % to preprocess for 1.5s after stimulus onset. 
        if SUBJ ==10 && B==3 || SUBJ==10 && B ==8; 
            trl(end,:) = []; 
        end
        
        %------------------------------------------------------------------
        % PREPROCESSING
        cfg = [];
        cfg.dataset     = fullname; % EEG data name
        cfg.trl = trl;
        cfg.demean     = 'no';
        cfg.detrend = 'no';
        cfg.polyremoval   = 'no';
        cfg.lpfilter   = 'yes';                                                 % apply lowpass filter
        cfg.lpfreq     = 90;                                                    % lowpass at 80 Hz.
        cfg.hpfilter   = 'yes';                                                 % apply highpass filter
        cfg.hpfreq     = ARG.Highpass;
        cfg.hpfiltord = 4;
        cfg.reref         = 'no';                                               % referencing
        cfg.continuous = 'yes';
        data = ft_preprocessing(cfg);
        [dataX] = eegck_BiosemiEselect(data);                                   % just pick the electrodes and eye channels
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % RESAMPLE THE DATA
        cfg            = [];
        cfg.resamplefs = ARG.Resample;                                          % resampling freq
        cfg.detrend    = 'no';
        dataX           = ft_resampledata(cfg, dataX);
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % PLOT RAW DATA BEFORE ICA
        XXtemp = zeros(length(dataX.trial),length(dataX.time{1}));
        XXtemp = eegck_trials2mat(dataX); 
        imagesc(sq(mean(XXtemp,2)))
        
%         for ii = 1:length(dataX.trial)
%             XXtemp(ii,:) = mean(dataX.trial{ii}(:,:),1);                        % average over trials
%         end
%         figure(),plot(dataX.time{1},XXtemp');
%         xlabel('Time'); ylabel('Amp');
%         suptitle(sprintf('S0%d Before ICA, Block %d',SUBJ,B));
%         saveppt(sprintf('Prepro S0%d',SUBJ));
%         close
        %------------------------------------------------------------------
        
        
        % COMPUTE EYE MOVEMENT SIGNALS
        % (THIS ASSUMES THAT ALL FOUR EOG ELECTRODES WERE USED)
        selchan = ft_channelselection({'EXG1' 'EXG2' 'EXG3' 'EXG4'}, dataX.label);
        data_eog = ft_selectdata(dataX, 'channel', selchan);
        % VEOG as difference between (EX3 - EX1) and (EX4-EX2)
        % HEOG as difference EX3-EX4
        % REOG as difference mean(EX1-4) - mean (A18,A20 A31)
        
%         %         selchan = eegck_returnEIndex({'A18','A20','A31'},dataX.label);
%         selchan = eegck_returnEIndex({'C17','C21','B24'},dataX.label);
%         for t=1:length(data_eog.trial)
%             veog = (data_eog.trial{t}(3,:)- data_eog.trial{t}(1,:))+ ...
%                 (data_eog.trial{t}(4,:)- data_eog.trial{t}(2,:));
%             veog = veog/2;
%             heog = (data_eog.trial{t}(3,:)- data_eog.trial{t}(4,:));
%             reog = mean(data_eog.trial{t},1)-mean(dataX.trial{t}(selchan,:));
%             % z-score each
%             data_eog.trial{t} = [cksig2z(veog);cksig2z(heog);cksig2z(reog)];
%         end
%         data_eog.label ={'VEOG','HEOG','REOG'};
%         %------------------------------------------------------------------
        
        selchan = [20 18 31 17 20 30];
        for t=1:length(data_eog.trial)
            veog = (data_eog.trial{t}(3,:)- data_eog.trial{t}(1,:))+(data_eog.trial{t}(4,:)- data_eog.trial{t}(2,:));
            veog = veog/2;
            heog = (data_eog.trial{t}(3,:)- data_eog.trial{t}(4,:));
            reog = mean(data_eog.trial{t},1)-mean(data.trial{t}(selchan,:));
            % z-score each
            data_eog.trial{t} = [cksig2z(veog);cksig2z(heog);cksig2z(reog)];
        end
        data_eog.label ={'VEOG','HEOG','REOG'};
    
    
        
        %------------------------------------------------------------------
        % ICA ARTIFACT REJECTION
        cfg            = [];
        cfg.method = 'runica';
        switch ARG.Remove
            case {'fast'}                                                       % do PCA to speed up process if searching only for frontal topos
                cfg.runica.pca = 40;
        end
        cfg.runica.maxsteps = 130;
        cfg.channel = dataX.label(1:128);                                       % this makes sure only eeg channels are used
        comp           = ft_componentanalysis(cfg, dataX);
        % close all
        
%         % DISPLAY COMPONENTS
        figure('units','normalized','outerposition',[0 0 1 1])
        cfg = [];
        cfg.component = [1: length(comp.label)];                                % specify the component(s) that should be plotted
        cfg.layout    = 'biosemi128.lay';                                       % specify the layout file that should be used for plotting
        cfg.comment   = 'no';
        ft_topoplotIC(cfg, comp)
%         
        % COMPUTE SUPPORT OF EACH COMPONENT
        s3 = comp.topo;
        s3 = s3-repmat(min(s3,[],1),[128,1]);
        s3 = s3./repmat(max(s3,[],1),[128,1]);
        s3= sum(s3);
        ARG.ICA.support = s3;
        
        % AUTOMATIC COMPONENT SELECTION BASED ON TOPOGRAPHY
        Artif_temp = eegck_artiftemp_biosemi(128);
        art_corr = [];
        for l=1:size(comp.topo,2)
            for k=1:size(Artif_temp,2)
                art_corr(l,k) = corr(comp.topo(:,l),Artif_temp(:,k));
            end
        end
        remove_eye = find( max(art_corr(:,(1:3)),[],2)>ARG.ICA_artif_thr);
        remove_all = find(max(art_corr(:,(4:end)),[],2)>ARG.ICA_artif_thr);
        ARG.ICA.art_corr = art_corr;
        
        % POWER SPECTRUM COMPONENTS
        % compute power spectra of components and removme components with a low
        % ratio between low-frequency power and high frequency power
        switch ARG.Remove
            case {'all'}
                [Power,fax,ratio] = eegck_componentPSD(comp);
                remove_power = find(ratio<ARG.ICA_spectrum_thr);
                good = find(ratio>ARG.ICA_spectrum_thr);
                % check those with strange power and remove if
                % correlations with templates
                remove2 =  find(max(abs(art_corr(remove_power,:)),[],2)>0.5);
                remove_all = [remove_all;remove_power(remove2)];
                ARG.ICA.ratio = ratio;
        end
        
        % EYE MOVEMENT/COMPONENT CORRELATIONS
        % compute correlations with eye movements and suggest components absed
        % on this
        switch ARG.Remove
            case {'all'}
                EyeAnalysis =  eegck_analyseEOG(ft_selectdata(data_eog, 'channel',{'VEOG','HEOG','REOG'}),0);
                M = eegck_trials2mat(comp);
                EyeAnalysis.Saccade  = EyeAnalysis.Saccade';
                EyeAnalysis.Movement  = EyeAnalysis.Movement';
                CC=[];
                for e=1:size(M,1)
                    tmp = (sq(M(e,:,:)));
                    tmp = ck_filt(tmp,dataX.fsample,[30,80],'band',6);
                    tmp = abs(ck_filt_hilbert(tmp))';
                    CC(e,1) = corr(tmp(:),EyeAnalysis.Saccade(:));
                    CC(e,2) = corr(tmp(:),EyeAnalysis.Movement(:));
                end
                ARG.ICA.all_eyecorr = CC;
                thr = std(CC(:));
                CC = max(abs(CC),[],2);
                remove_eye_corr = find((CC>ARG.ICA_eyecorr));
                ARG.ICA.remove_eye_corr = remove_eye_corr;
                EyeAnalysis=[];
                M=[];
        end
        
        % SUGGESTED COMPONENTS TO REMOVE
        switch ARG.Remove
            case {'fast'}
                % remove frontal topos only - automatic process
                remove = remove_eye;
                ARG.ICA.not_removed=[];
            case {'all'}
                % remove based on topo
                remove = [remove_eye;remove_all;remove_eye_corr];
                % find those with seem to have strange power but are not rejected and display
                ARG.ICA.not_removed  = setdiff(remove_power,remove);
        end
        remove = unique(remove);
        remove = sort(remove);
        fprintf('Suggested components to remove:  \n');
        for l=1:length(remove)
            fprintf('%d \n ',remove(l));
        end
        fprintf('\n');
        fprintf('%d \n', size(remove,1));
        %------------------------------------------------------------------
        
        
%         %------------------------------------------------------------------
%         % DISPLAY COMPONENTS
%         figure(1);
%         chld = get(gcf,'Children');
%         ncomp = length(comp.label);
%         for l=1:length(remove)
%             set(get(chld(ncomp-remove(l)+1),'Title'),'Color',[1 0 0 ])
%         end
%         % red = removed
%         for l=1:length(ARG.ICA.not_removed)
%             set(get(chld(ncomp-ARG.ICA.not_removed(l)+1),'Title'),'Color',[0 1 0 ])
%         end
% %         saveppt(sprintf('Prepro S0%d',SUBJ));
% %         close
%         %------------------------------------------------------------------
        
        
%         
%         %------------------------------------------------------------------
%         % DISPLAY SPECTRA AND COMPONENTS
%         % display spectra and components that may be conspicous but were not removed
%         switch ARG.Remove
%             case {'all'}
%                 figure(4);
%                 set(gcf,'Position',[ 0 0        1484        1000]);
%                 n = numel( ARG.ICA.not_removed)*2;
%                 nyplot = ceil(sqrt(n));  nxplot = ceil(n./nyplot);
%                 load('FtDummy');
%                 eegck_topcfg;
%                 for k=1:length( ARG.ICA.not_removed)
%                     subplot(nxplot, nyplot,k*2-1); hold on
%                     plot(fax,(log(Power(good,:))),'k');
%                     plot(fax,log(Power( ARG.ICA.not_removed(k),:)),'r','LineWidth',2);
%                     axis([fax(1) fax(end) min(log(Power(:))) max(log(Power(:)))]);
%                     title(sprintf('component %d', ARG.ICA.not_removed(k)));
%                     subplot(nxplot, nyplot,k*2);
%                     EvpDummy.avg = comp.topo(:, ARG.ICA.not_removed(k));
%                     ft_topoplotER(cfg_topo,EvpDummy);
%                 end
% %                 suptitle(sprintf('S0%d,Block %d',SUBJ,B));
% %                 saveppt(sprintf('Prepro S0%d',SUBJ));
%         end
%         close all
%         %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % REMOVE COMPONENTS
        cfg = [];
        cfg.component = remove;
        dataX = ft_rejectcomponent(cfg, comp);
        ARG.ICA.removed_components = remove;
        % just preserve the data before removing trials;
        dataTemp = dataX;
        %------------------------------------------------------------------
%         
%         % check for bad channels
%         out = eegck_maxampft(dataX,(1:128));
%         imagesc(out);
%         addpath('W:\CIRCLE_EXP\LOG\EEG');
        
        
% % %         %------------------------------------------------------------------
% % %         % REPAIR ANY BAD CHANNELS IN SUBJECTS
% % % % %         fprintf('\n\n')
% % % %        
% dataX = dataTemp; 
% % 
% channels = [25 26 27 7 72 73 27 28 39 40 92 93 94]; 
% badchannels = dataX.label(channels);
% cfg = [];
% cfg.method = 'nearest';
% cfg.badchannel = badchannels;
% cfg.layout = 'biosemi128.lay';
% cfg.trials  = 'all';
% cfg.sens = 'biosemi128.lay';
% cfg_neighb = eegsb_BiosemiNbrs128();
% cfg.neighbours = cfg_neighb.neighbours;
% [interp] = ft_channelrepair(cfg, dataX);
% dataX = interp; % interpolated data
% fprintf('Removed Bad Channels \n');
% %------------------------------------------------------------------
% % % % %
%                %------------------------------------------------------------------
        % PLOT RAW DATA BEFORE ICA
        XXtemp = zeros(length(dataX.trial),length(dataX.time{1}));
        XXtemp = eegck_trials2mat(dataX); 
        imagesc(sq(mean(XXtemp,2)),[-3 3]) 
% % % % %         

% %         
        %------------------------------------------------------------------
        % REMOVE TRIALS WITH EXCESSIVE AMPLITUDE
        out = eegck_maxampft(dataX,(1:128));
        cfg = [];
        cfg.trials = find(sum(abs(out)>120,1)==0)
%         
        %%
        dataX = ft_redefinetrial(cfg, dataX);
        fprintf('\n\n')
        data_eog = ft_redefinetrial(cfg,data_eog);
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % APPEND EOG AND EEG DATA
        dataX = ft_appenddata([], dataX, data_eog);
        %------------------------------------------------------------------
        
        
%         ------------------------------------------------------------------
        % REPLOT AFTER ICA
        XXtemp = zeros(length(dataX.trial),length(dataX.time{1}));
        XXtemp = eegck_trials2mat(dataX); 
%         
%         for ii = 1:length(dataX.trial)
%             XXtemp(ii,:) = mean(dataX.trial{ii}(:,:),1);
%         end
%         figure(),plot(dataX.time{1},XXtemp');
%         xlabel('Time'); ylabel('Amp');
%         box off
%         title(sprintf('S0%d After ICA, Block %d',SUBJ,B));
%         saveppt(sprintf('Prepro S0%d',SUBJ));
%         close all
%         %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % save data
        sname = sprintf('PreproICA_S0%d_B%d.mat',SUBJ,B);
        save(sname,'dataX','A','ARG');
        fprintf('\n\n')
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        % EOG ANALYSIS
        EyeAnalysis =  eegck_analyseEOG(ft_selectdata(dataX, 'channel',{'VEOG','HEOG','REOG'}),0);
        j = find( (dataX.time{1}>=-0.5).*(dataX.time{1}<=1.2));             % select time region of interest – adapdt your paradigm e,g, -0.5 to 1.2
        MaskM =EyeAnalysis.Movement(:,j);                                   % sum over time
        MaskS =EyeAnalysis.Saccade(:,j);
        EyeAnalysis=[];
        
        % two thresholds to select trials based on number of bad time points
        ARG.EOGS = 11;
        ARG.EOGM = 20;
        
        % clean signal
        cfg = [];
        cfg.trials = find((sum(MaskM,2)<ARG.EOGM).*(sum(MaskS,2)<ARG.EOGS))
        %%
        dataX = ft_redefinetrial(cfg, dataX);
        ARG.EOGremove = cfg;
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        newName = sprintf('PreproICA_S0%d_B%d_EOG_out.mat',SUBJ,B);
        save(newName,'dataX','ARG','A');
        % fprintf('\n\n\n\');
        fprintf('Done:  Preprocessing S0%d Block %d  \n',SUBJ,B)
        %------------------------------------------------------------------
        
        %%
    end
    
    
    
    %------------------------------------------------------------------
    % CONCATENATE BLOCKS
    files = dir('PreproICA_*_EOG_out.mat');
    
    dataNew = [];
    for k = 1:length(files);
        load(files(k).name,'dataX');
        dataNew{k} = dataX;
    end
    dataX = ft_appenddata([],dataNew{:});
    
    save(sprintf('Prepro_all_S0%d_EOG_out',SUBJ),'dataX')
    %------------------------------------------------------------------
    
    
    %%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY NOISE CLEANING FOR EACH BLOCK - INDIVIDUAL PARTICIPANTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dataList = dir('PreproICA_S*_EOG_out.mat');
    c=1;co=1;
    nfiles = length(dataList);
    
    % Get data from the subject from the preprocessing file
    for file=1:nfiles
        sname = dataList(file).name;
        Data = load(sname,'ARG');
        Stats_sup(co,:) = Data.ARG.ICA.support;
        Stats_sup2(co) = median(Data.ARG.ICA.support(Data.ARG.ICA.removed_components));
        Stats_thr(co,:,:) = Data.ARG.ICA.art_corr;
        Stats_ratio(co,:) = Data.ARG.ICA.ratio;
        Stats_eyec(co,:,:) = Data.ARG.ICA.all_eyecorr;
        Stats_remove(co,:) = [length(Data.ARG.ICA.remove_eye_corr),length(Data.ARG.ICA.removed_components)];
        Sub_list(co,:) = [SUBJ,file];
        co=co+1;
    end
    
    cvec = 'rgbkmcrgbkmcrgbkmcrgbkmcrgbkmcrgbkmc';
    
    %--------------------------------------------------------------------------
    % PLOT RESULTS
    N = size(dataList,1);
    L =  (1 :length(dataList)); % axis length
    figure('units','normalized','outerposition',[0 0 0.7 0.8])
    
    % correlation with templates
    subplot(2,2,1); hold on;
    for k=1:N
        tmp = max(abs(sq(Stats_thr(k,:,:))),[],2);
        plot(k,tmp,[cvec(Sub_list(k,2)) 'o'])
    end
    title('Template correlation');
    xlim([0 length(Sub_list)+1])
    set(gca,'XTick',L,'XTickLabel',Sub_list(L,2),'TickDir','out')
    line([0 N+1],[Data.ARG.ICA_artif_thr Data.ARG.ICA_artif_thr])
    xlabel('block')
    
    % Power ratio
    subplot(2,2,2); hold on;
    for k=1:N
        tmp =Stats_ratio(k,:);
        plot(k,tmp,[cvec(Sub_list(k,2)) 'o'])
    end
    title('Power ratio')
    axis([0 N+1 0 30]);                                                 % we care most abotu those with low ratio
    line([0 N+1],[Data.ARG.ICA_spectrum_thr Data.ARG.ICA_spectrum_thr])
    set(gca,'XTick',L,'XTickLabel',Sub_list(L,2),'TickDir','out')
    xlabel('block')
    
    % Eye Correlation
    subplot(2,2,3); hold on;
    for k=1:N
        tmp =Stats_eyec(k,:);
        plot(k,tmp,[cvec(Sub_list(k,2)) 'o'])
    end
    title('Eye movement correlation')
    xlim([0 length(Sub_list)+1])
    line([0 N+1],[Data.ARG.ICA_eyecorr Data.ARG.ICA_eyecorr])
    set(gca,'XTick',L,'XTickLabel',Sub_list(L,2),'TickDir','out')
    xlabel('block')
    
    % Number of removed components
    subplot(2,2,4); hold on;
    for k=1:N
        plot(k,Stats_remove(k,1),[cvec(Sub_list(k,1)) 'o'])
        plot(k,Stats_remove(k,2),[cvec(Sub_list(k,1)) 's'])
    end
    plot(Stats_remove(:,1),'k');
    plot(Stats_remove(:,2),'r');
    xlim([0 length(Sub_list)+1])
    title('# removed components (total/eyeC)');
    set(gca,'XTick',L,'XTickLabel',Sub_list(L,2),'TickDir','out')
    xlabel('block');
    
    suptitle(sprintf('S0%d',SUBJ));
    saveppt(sprintf('Prepro S0%d',SUBJ));
    close all
    clc
    fprintf('S0%d Done .... \n',SUBJ)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end

% cd('W:\CIRCLE_EXP\LOG\EEG');
