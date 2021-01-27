%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS SCRIPT: STEPH B PHD PROJECT 2 (8/2/17)
% AV CROSSMODAL CORRESPONDENCES: High/low tone or small/large circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEHAVIOURAL ANALYSIS: INDIVIDUAL SUBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
cd('W:/CIRCLE_EXP/LOG/EEG')
subs = dir('Subj_*');
CUTOFF = [0.3 1.2];                                                                 % reaction time cutoff

for SUBJ = 1:length(subs);
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',SUBJ))
    
    % different trials
    load(sprintf('S0%d_allT.mat',SUBJ))
    trialinfo = allT;
    trialinfo = trialinfo(trialinfo(:,7)>CUTOFF(1) & trialinfo(:,7)<CUTOFF(2),:);
    clear AC AIC VC VIC AClog AIClog VClog VIClog
    
    % reaction times
    RT = trialinfo(:,7);
    logRT = log(1./RT);
    
    % indices
    ACind = find(trialinfo(:,2)<=2 & trialinfo(:,3)==1,7);
    AICind = find(trialinfo(:,2)<=2 & trialinfo(:,3)==2,7);
    VCind = find(trialinfo(:,2)>=3 & trialinfo(:,3)==1,7);
    VICind = find(trialinfo(:,2)>=3 & trialinfo(:,3)==2,7);
    
    % normal RT for different congruency conditions
    AC = trialinfo(trialinfo(:,2)<=2 & trialinfo(:,3)==1,7);
    AIC = trialinfo(trialinfo(:,2)<=2 & trialinfo(:,3)==2,7);
    VC = trialinfo(trialinfo(:,2)>=3 & trialinfo(:,3)==1,7);
    VIC = trialinfo(trialinfo(:,2)>=3 & trialinfo(:,3)==2,7);
    
    % Log transformed RT
    AClog = log(1./AC);
    AIClog = log(1./AIC);
    VClog = log(1./VC);
    VIClog = log(1./VIC);
    
    % store the different reaction times in variables
    RTcond = {AC AIC VC VIC};                                                       % normal RT
    RTlog = {AClog AIClog VClog VIClog};                                            % log RT
    
    % work out means of the different RT analyses RT
    RTG = [mean(AC) mean(AIC) mean(VC) mean(VIC)];
    RTG_log = [mean(AClog) mean(AIClog) mean(VClog) mean(VIClog)];
    
    % save stuff
    save(sprintf('S0%d_BEH_CUTOFF.mat',SUBJ),'trialinfo',...
        'RT','logRT','RTcond','RTlog','AC','AIC','VC','VIC',...
        'ACind','AICind','VCind','VICind')
    
    fprintf('S0%d Done \n',SUBJ);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take out trials that were cut out of behaviorual files for RT from the
% EEG files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
CUTOFF = [0.3 1.2];             % reaction time cutoff

for SUBJ = 1:length(subs);
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',SUBJ))
    load(sprintf('Prepro_all_S0%d_EOG_out.mat',SUBJ))
    
    % save old data before cutoff
    save(sprintf('S0%d_EEG_ALL.mat',SUBJ),'dataX')
    
    % new data (trials with extreme reaction times taken out
    cfg = [];
    cfg.trials =  find(dataX.trialinfo(:,7)>CUTOFF(1) & dataX.trialinfo(:,7)<CUTOFF(2));
    dataX = ft_redefinetrial(cfg,dataX);
    
    save(sprintf('S0%d_EEG_CUTOFF.mat',SUBJ),'dataX')
    fprintf('S0%d Done \n',SUBJ)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECODING ANALYSIS: INDIVIDUAL SUBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
cd('W:/CIRCLE_EXP/LOG/EEG');
subs = dir('Subj_*');
addpath('W:\CIRCLE_EXP\LOG\EEG')
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')

% DECODER SETTINGS
ARG.AmpThr = 100;                                                                   % amplitude threshold cutoff
ARG.Twin = [-0.3 1.0];                                                              % time points for decoder
ARG.Tsteps = 5;                                                                     % every nth time point  (5ms)
ARG.winL = 7;                                                                       % decoder length        (7 = 30ms)

% decoding parameters
cfg_decode =[];
cfg_decode.CVsteps = 10;                                                            % cross-validation steps
cfg_de.channels = (1:128);                                                          % channels
cfg_decode.reg.gamma = 0.1;                                                         % regularisation parameter

% Decoder prepro settings
cfg_prepro = [];
cfg_prepro.lpfilter   = 'yes';
cfg_prepro.lpfreq     = 30;
cfg_prepro.hpfilter   = 'yes';
cfg_prepro.hpfreq     = 1;
cfg_prepro.demean ='no';
cfg_prepro.reref = 'no';

[labels,index,labelGroup] = eegsb_GroupLabels128();                                 % Group the labels into 4
load('EVP128.mat')                                                                  % Dummy EEG file to use fieldtrip structure

%%
% LOOP THROUGH SUBJECTS
for SUBJ = 2:length(subs);
    
    close all
    % For the Topoplots: dummy file
    clearvars -except CUTOFF subs ARG cfg_decode cfg_prepro labels index labelGroup EvpDummy SUBJ
    
    ARG.Twin = [-0.3 1.0];          % time points for decoder
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',SUBJ))
    load(sprintf('S0%d_EEG_CUTOFF.mat',SUBJ),'dataX')
    timings = dataX.time{1};
    trange = find( (timings>=ARG.Twin(1)).*(timings<=ARG.Twin(end)));
    Decode_times = timings(trange);
    
    % all trials
    dataTarget = ft_preprocessing(cfg_prepro,dataX);
    MT = eegck_trials2mat(dataTarget);
    trialinfo = dataTarget.trialinfo;
    
%     %----------------------------------------------------------------------
%     % DECODER 1: AUDITORY FROM VISUAL
%     %----------------------------------------------------------------------
%     YCon = trialinfo(:,2)<=2;                                                       % 1 = AUD, 0 = VIS
%     Az = zeros(1,length(trange));
%     A = zeros(128,length(trange));
%     Y = zeros(length(YCon),length(trange));
%     Wsave = zeros(length(trange),128);
%     Csave = zeros(1,length(trange));
%     YT = zeros(size(MT,2),length(trange));
%     Az_YT = zeros(1,length(trange));
%     
%     % loop through each time point in the trial
%     for t=1:length(trange)
%         act2 = mean(MT(1:128,:,trange(t)+(1:ARG.winL)),3);
%         [A(:,t),Az(t),Y(:,t),Wmix]= eegck_LDA(cfg_decode,mean(MT(1:128,:,trange(t)+(1:ARG.winL)),3)',YCon);
%         Cmix = Wmix(end); Wmix = Wmix(1:end-1);
%         Wsave(t,:) = Wmix; Csave(t) = Cmix;
%     end
%     fprintf('S0%d SENSORY DECODER DONE... \n',SUBJ)
%     
%     
%     % PLOTS: SENSORY DECODER
%     figure('units','normalized','outerposition',[0 0 1 1])
%     stimWind = find(Decode_times>=-0.1 & Decode_times<=0.6);                        % just plot a smaller window rather than whole trial
%     
%     % Plot Az value
%     subplot 221
%     plot(Decode_times,Az)
%     vline(0,':k'); vline(0.3,':k'); hline(0.5,':k')
%     title('AUD VS VIS');
%     
%     % Topoplot: use the scalp projection from the highest Az value
%     subplot 222
%     cfg = [];
%     cfg.layout = 'biosemi128.lay';
%     cfg.comment = ' ';
%     cfg.commentpos = 'title';
%     cfg.colorbar = 'yes';
%     [~,J2] = max(Az);
%     EvpDummy.avg = sq(A(:,J2));
%     ft_topoplotER(cfg,EvpDummy)
%     title(sprintf('AUD vs VIS %s',num2str(Decode_times(J2))))
%     
%     % Y Decoding signal plots
%     aud = Y(YCon==1,:);                                                             % AUDITORY TRIALS
%     vis = Y(YCon==0,:);                                                             % VISUAL TRIALS
%     AC = Y(trialinfo(:,2)<=2 & trialinfo(:,3)==1,:);                                % AUDITORY CONGRUENT TRIALS
%     AIC = Y(trialinfo(:,2)<=2 & trialinfo(:,3)==2,:);                               % AUDITORY INCONGRNT TRIALS
%     VC = Y(trialinfo(:,2)>=3 & trialinfo(:,3)==1,:);                                % VISUAL CONGRUENT TRIALS
%     VIC = Y(trialinfo(:,2)>=3 & trialinfo(:,3)==2,:);                               % VISUAL INCONGRNT TRIALS
%     
%     % Plot Y signals, not split up by congruency:
%     subplot 223; hold on
%     plot(Decode_times,sq(mean(aud)),'LineWidth',2);
%     plot(Decode_times,sq(mean(vis)),'m','LineWidth',2);
%     xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%     vline(0,':k'); vline(0.3,':k'); hline(0.5,':k')
%     title('AUD VS VIS');
%     
%     % Plot Y signals, split up by congruency:
%     subplot 224; hold on
%     plot(Decode_times,sq(mean(AC)),'LineWidth',2);
%     plot(Decode_times,sq(mean(AIC)),'k','LineWidth',2);
%     plot(Decode_times,sq(mean(VC)),'r','LineWidth',2);
%     plot(Decode_times,sq(mean(VIC)),'m','LineWidth',2);
%     legend({'AC','AIC','VC','VIC'},'FontSize',6)
%     xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%     vline(0,':k'); vline(0.3,':k'); hline(0.5,':k')
%     title('AUD VS VIS');
%     
%     suptitle(sprintf('S0%d SENSORY (A vs V) DECODER \n %s',SUBJ))
%     legend({'AC','AIC','VC','VIC'},'FontSize',6)
%     saveppt(sprintf('S0%d',SUBJ))
%     
    
    %----------------------------------------------------------------------
    % DECODER 2: AUDITORY HIGH VS LOW, VISUAL SMALL VS LARGE
    %----------------------------------------------------------------------
    condition = trialinfo(:,2);                                                     % condition 1:4 (1 = AUD H, 2 = AUD L
    aud = find(condition<=2);                                                       % auditory trials
    vis = find(condition>=3);                                                       % visual trials
    trialinfo(aud,1)=1;                                                             % add in whether a trial is auditory or visual
    trialinfo(vis,1)=2;
    
    % Stimulus index
    STIM{1} = find(trialinfo(:,2)==1);                                              % high tone
    STIM{2} = find(trialinfo(:,2)==2);                                              % low tone
    STIM{3} = find(trialinfo(:,2)==3);                                              % small circle
    STIM{4} = find(trialinfo(:,2)==4);                                              % large circle
    
    % AUDITORY TRIALS - HIGH VS LOW
    YConA = ismember(aud,STIM{1});                                                  % 1 = high tone, 0 = low tone.
    Az_aud = zeros(1,length(trange));
    A_aud = zeros(128,length(trange));
    Yaud = zeros(length(YConA),length(trange));
    YTA = zeros(size(MT,2),length(trange));
    TopoAUD = zeros(128,length(trange));
    WsaveA = zeros(length(trange),128);
    
    % VISUAL TRIALS - HIGH VS LOW
    YConV = ismember(vis,STIM{3});                                                  % 1 = small circle, 0 = large circle
    Az_vis = zeros(1,length(trange));
    A_vis = zeros(128,length(trange));
    Yvis = zeros(length(YConV),length(trange));
    YTV = zeros(size(MT,2),length(trange));
    TopoVIS = zeros(128,length(trange));
    WsaveV = zeros(length(trange),128);
    CsaveA = zeros(1,length(trange));
    CsaveV = zeros(1,length(trange));
    
    
    % USE CURRENT TIME POINT WEIGHT
    for t=1:length(trange)
        
        % Decode auditory
        [A_aud(:,t),Az_aud(t),Yaud(:,t),Wmix]= eegck_LDA(cfg_decode,mean(MT(1:128,aud,trange(t)+(1:ARG.winL)),3)',YConA);
        Cmix = Wmix(end); Wmix = Wmix(1:end-1);
        WsaveA(t,:) = Wmix; CsaveA(t) = Cmix;
        
        % Decode visual
        [A_vis(:,t),Az_vis(t),Yvis(:,t),Wmix]= eegck_LDA(cfg_decode,mean(MT(1:128,vis,trange(t)+(1:ARG.winL)),3)',YConV);
        Cmix = Wmix(end);
        Wmix = Wmix(1:end-1);
        WsaveV(t,:) = Wmix;
        CsaveV(t) = Cmix;
        
    end
    fprintf('S0%d Normal decoder done... \n',SUBJ);
    
    
%     % Save sensory decoder
%     SD.Az = Az;
%     SD.A = A;
%     SD.Y = Y;
%     SD.YCon = YCon;
%     SD.Decode_times = Decode_times;
    
    % Save non projected decoder
    ND.AudAz = Az_aud;
    ND.TopoAud = A_aud;
    ND.Yaud = Yaud;
    ND.VisAz = Az_vis;
    ND.TopoVis = A_vis;
    ND.Yvis = Yvis;
    ND.YConA = YConA;
    ND.YConV = YConV;
    ND.Decode_times = Decode_times;
    ND.WsaveA = WsaveA;
    ND.WsaveV = WsaveV;
    
%     save(sprintf('S0%d_DECODE_CUTOFF_2.mat',SUBJ),'SD','ND','Decode_times','trange',...
%         'trialinfo','YConA','YConV','YCon')
        save(sprintf('S0%d_DECODE_CUTOFF_2.mat',SUBJ),'ND','Decode_times','trange',...
        'trialinfo','YConA','YConV')
%     
%     
%     
%     %----------------------------------------------------------------------
%     % Plots: HIGH VS LOW TONE - NON PROJECTED TRIALS
%     clear Y
%     data = ND;
%     name = {'AUD H vs. L','VIS S vs L'};
%     sname = {'AUD H vs L, VIS S vs L, NON PROJECTED trials'};
%     clear AHC AHIC ALC ALIC VSC VSIC VLC VLIC
%     
%     Az(1,:) = data.AudAz;
%     Az(2,:) = data.VisAz;
%     Topo(:,:,1) = data.TopoAud;
%     Topo(:,:,2) = data.TopoVis;
%     Y{1} = data.Yaud;
%     Y{2} = data.Yvis;
%     trials = trialinfo;
%     
%     AUD = find(trials(:,2)<=2);                                                     % AUDITORY
%     VIS = find(trials(:,2)>=3);                                                     % VISUAL
%     AHC = Y{1}(trials(AUD,2)==1 & trials(AUD,3)==1,:);                              % AUDITORY HIGH, CONGRUENT
%     AHIC = Y{1}(trials(AUD,2)==1 & trials(AUD,3)==2,:);                             % AUDITORY HIGH, INCONGRNT
%     ALC = Y{1}(trials(AUD,2)==2 & trials(AUD,3)==1,:);                              % AUDITORY LOW, CONGRUENT
%     ALIC = Y{1}(trials(AUD,2)==2 & trials(AUD,3)==2,:);                             % AUDITORY LOW, INCONGRNT
%     
%     VSC = Y{2}(trials(VIS,2)==3 & trials(VIS,3)==1,:);                              % VISUAL SMALL, CONGRUENT
%     VSIC = Y{2}(trials(VIS,2)==3 & trials(VIS,3)==2,:);                             % VISUAL SMALL, INCONGRNT
%     VLC = Y{2}(trials(VIS,2)==4 & trials(VIS,3)==1,:);                              % VISUAL LARGE, CONGRUENT
%     VLIC = Y{2}(trials(VIS,2)==4 & trials(VIS,3)==2,:);                             % VISUAL LARGE, INCONGRNT
%     cc = distinguishable_colors(5);
%     
%     % Plot Az trials
%     for k = 1:2;
%         subplot 231; hold on
%         plot(Decode_times,Az(k,:),'Color',cc(k,:),'LineWidth',2);
%         xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%         ylim([0.4 1])
%         vline(0,':k'); vline(0.3,':k'); hline(0.5,':k')
%     end
%     legend(name)
%     
%     % Plot Topography
%     np = 2; load('EVP128.mat')
%     cfg = [];
%     cfg.layout = 'biosemi128.lay';
%     cfg.comment = ' ';
%     cfg.commentpos = 'title';
%     cfg.colorbar = 'yes';
%     for k = 1:2;
%         subplot(2,3,np)
%         [~,J2] = max(Az(k,stimWind));
%         EvpDummy.avg = sq(Topo(:,J2,k));
%         ft_topoplotER(cfg,EvpDummy)
%         title(sprintf('%s %s',name{k},num2str(Decode_times(stimWind(J2)))))
%         np = np+1;
%     end
%     
%     % Y SIGNALS AUD
%     subplot 234; hold on
%     plot(Decode_times,sq(mean(AHC)),'Color','b','LineWidth',2);                     % CONGRUENT
%     plot(Decode_times,sq(mean(AHIC)),'Color',[0.5 0.85 1],'LineWidth',2)            % INCONGRUENT
%     plot(Decode_times,sq(mean(ALC)),'Color',[0.7 0.2 1],'LineWidth',2)              % CONGRUENT
%     plot(Decode_times,sq(mean(ALIC)),'Color',[0.7 0.7 1],'LineWidth',2)             % INCONGRUENT
%     xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%     vline(0,':k'); vline(0.3,':k'); hline(0,':k')
%     legend({'AHC','AH IC','ALC','AL IC'},'FontSize',6)
%     title('AUD')
%     
%     % Y SIGNALS VIS
%     subplot 235; hold on
%     plot(Decode_times,sq(mean(VSC)),'Color',[1 0.1 0.4],'LineWidth',2)              % CONGRUENT
%     plot(Decode_times,sq(mean(VSIC)),'Color',[1 0.78 0.80],'LineWidth',2)           % INCONGRUENT
%     plot(Decode_times,sq(mean(VLC)),'Color',[1 0.6 0.1],'LineWidth',2)              % CONGRUENT
%     plot(Decode_times,sq(mean(VLIC)),'Color',[1 0.85 0.4],'LineWidth',2)            % INCONGRUENT
%     xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%     vline(0,':k'); vline(0.3,':k'); hline(0,':k')
%     legend({'VSC','VS IC','VLC','VL IC'},'FontSize',6)
%     title('VIS')
%     
%     % DIFFERENCE PLOT
%     ADIFF(1,:) = mean(AHC) - mean(ALC);
%     ADIFF(2,:) = mean(AHIC) - mean(ALIC);
%     VDIFF(1,:) = mean(VLC) - mean(VSC);
%     VDIFF(2,:) = mean(VLIC) - mean(VSIC);
%     subplot 236; hold on;
%     plot(Decode_times,ADIFF(1,:),'Color','b','LineWidth',2)                         % AUD CONG
%     plot(Decode_times,ADIFF(2,:),'Color',[0.5 0.85 1],'LineWidth',2)                % AUD INCONG
%     plot(Decode_times,VDIFF(1,:)*-1,'Color',[1 0.1 0.4],'LineWidth',2)
%     plot(Decode_times,VDIFF(2,:)*-1,'Color',[1 0.78 0.80],'LineWidth',2)
%     xlim([Decode_times(stimWind(1)) Decode_times(stimWind(end))]);
%     hline(0,':k');vline(0,':k'); vline(0.3,':k');
%     
%     title('DIFF (AH-AL) and (VL-SC)')
%     
%     suptitle(sname)
%     legend({'AC','AIC','VC','VIC'},'FontSize',6)
%     saveppt(sprintf('S0%d',SUBJ))
%     %----------------------------------------------------------------------
%     
    
    fprintf('S0%d Decoders Done... \n',SUBJ)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REGRESSION ANALYSIS : INDIVIDUAL SUBJECTS, NOT SPLIT UP BY CONGRUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
cd('W:/CIRCLE_EXP/LOG/EEG');
subs = dir('Subj_*');
addpath('W:\CIRCLE_EXP\LOG\EEG')
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')
aw = zeros(length(subs),261);


for SUBJ = 1:length(subs);
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',SUBJ))
    load(sprintf('S0%d_BEH_CUTOFF.mat',SUBJ),'RT');
    load(sprintf('S0%d_DECODE_CUTOFF.mat',SUBJ));
    
    % NORMAL CURRENT TIME POINT DECODER
    TIMES = ND.Decode_times;
    audT = find(trialinfo(:,2)<=2);                                                 % auditory trials
    
    %----------------------------------------------------------------------
    % Auditory Trials first
    %----------------------------------------------------------------------
    indA = find(ND.YConA==0);                                                       % flip decoding signals that were assigned 0
    Y_AUD = ND.Yaud;
    YflipAUD = ND.Yaud;
    YflipAUD(indA,:) = YflipAUD(indA,:)*-1;
    RTuseA = RT(audT);                                                              % reaction times
    
    % build auditory regression model
    actA = YflipAUD;
    xA = zeros(size(actA,1),2);
    xA(1:length(RTuseA),1) = RTuseA;
    xA(:,2) = 1;
    
    % Regress reaction time against Y signal at each time point : all auditory trials
    AUD_Weight = zeros(length(TIMES),2);
    AUD_CI = zeros(length(TIMES),2);
    AUD_STATS = zeros(length(TIMES),4);
    for t = 1:length(ND.Decode_times);
        [b,CI,~,~,STATS] = regress(actA(:,t),xA);
        AUD_Weight(t,:) = b;                                                        % [coefficient, constant]
        AUD_CI(t,:) = CI(1,:);                                                      % confidence interval from coefficient
        AUD_STATS(t,:) = STATS;                                                     % [R^2 stat, F stat, pvalue, error variance estimate]
    end
    
    save(sprintf('S0%d_REG_NEW.mat',SUBJ),...
        'YflipAUD','Y_AUD','RTuseA','audT','AUD_Weight','AUD_CI','AUD_STATS',...
        'AC_ind','AIC_ind','AC','AIC');
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % Visual Trials second
    %----------------------------------------------------------------------
    visT = find(trialinfo(:,2)>=3);                                                 % visual trials
    indV = find(ND.YConV==0);
    Y_VIS = ND.Yvis;
    YflipVIS = ND.Yvis;
    YflipVIS(indV,:) = YflipVIS(indV,:)*-1;
    RTuseV = RT(visT);                                                              % reaction times
    
    % build visual regression model
    actV = YflipVIS;
    xV = zeros(size(actV,1),2);
    xV(1:length(RTuseV),1) = RTuseV;
    xV(:,2) = 1;
    
    % Regress reaction time against Y signal at each time point
    VIS_Weight = zeros(length(TIMES),2);
    VIS_CI = zeros(length(TIMES),2);
    VIS_STATS = zeros(length(TIMES),4);
    for t = 1:length(ND.Decode_times);
        [b,CI,~,~,STATS] = regress(actV(:,t),xV);
        VIS_Weight(t,:) = b;                                                        % [coefficient, constant]
        VIS_CI(t,:) = CI(1,:);                                                      % confidence interval from coefficient
        VIS_STATS(t,:) = STATS;                                                     % [R^2 stat, F stat, pvalue, error variance estimate]
    end
    
    save(sprintf('S0%d_REG_VIS.mat',SUBJ),...
        'YflipVIS','Y_VIS','RTuseV','visT','VIS_Weight','VIS_CI','VIS_STATS',...
        'VC_ind','VIC_ind','VC','VIC')
    
    fprintf('S0%d Regression Done \n',SUBJ);
    %----------------------------------------------------------------------]
    
end

cd('W:\CIRCLE_EXP\LOG\EEG');

%--------------------------------------------------------------------------
% Regression models above:
% 1 model | where Y activity auditory is regressed against Reaction time.
% Includes 1 preditor (RT) and constant (1).
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              GROUP ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% One subject (S19)'s data had to be discarded due to very noisy channels.
% They're left in the analysis above, just because that's how it started.
% But the group files below reads in subjects 1:18,20.

% Also, after the behavioural group file, the auditory and visual EEG
% trials are done separately in the code because I originally did them
% separately, and honestly can't be bothered integrating the two codes
% right now. They both work so it's fine.
% Auditory gets done first, then visual.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEHAVIOURAL ANALYSIS : GROUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG');
subs = [1:18 20];
n = length(subs);
Group_Perf = zeros(n,4);
Group_RT = zeros(n,4);
Group_RTL = zeros(n,4);

trialsAll = [];

for SUBJ = 1:n
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    load(sprintf('S0%d_BEH_CUTOFF',subs(SUBJ)),'RTcond','RTlog','trialinfo')
    
    %----------------------------------------------------------------------
    % Group Reaction Times
    %----------------------------------------------------------------------
    for k = 1:4;
        Group_RT(SUBJ,k) = mean(RTcond{k});                                         % normal
        Group_RTL(SUBJ,k) = mean(RTlog{k});                                         % log transformed
    end
    
    %----------------------------------------------------------------------
    % Group Performance Scores
    %----------------------------------------------------------------------
    % auditory congruent
    Group_Perf(SUBJ,1) = size(find(trialinfo(:,2)<=2 & trialinfo(:,3)==1 & ...
        trialinfo(:,8)==1),1)/size(find(trialinfo(:,2)<=2 & trialinfo(:,3)==1),1);
    
    % auditory incongruent
    Group_Perf(SUBJ,2) = size(find(trialinfo(:,2)<=2 & trialinfo(:,3)==2 & ...
        trialinfo(:,8)==1),1)/size(find(trialinfo(:,2)<=2 & trialinfo(:,3)==2),1);
    
    % visual congruent
    Group_Perf(SUBJ,3) = size(find(trialinfo(:,2)>=3 & trialinfo(:,3)==1 & ...
        trialinfo(:,8)==1),1)/size(find(trialinfo(:,2)>=3 & trialinfo(:,3)==1),1);
    
    % visual incongruent
    Group_Perf(SUBJ,4) = size(find(trialinfo(:,2)>=3 & trialinfo(:,3)==2 & ...
        trialinfo(:,8)==1),1)/size(find(trialinfo(:,2)>=3 & trialinfo(:,3)==2),1);
    
    fprintf('S0%d done \n',subs(SUBJ))
    trials{SUBJ} = trialinfo;
    trialsAll = cat(1,trialsAll,trialinfo);
end

cd('W:\CIRCLE_EXP\LOG\EEG\Files_Behaviour');
save('GROUP_BEH_CUTOFF_19.mat','Group_Perf','Group_RT','Group_RTL','trials','trialsAll');


%--------------------------------------------------------------------------
% Reaction Time Plots
%--------------------------------------------------------------------------
cc = [0.8 0.8 0.8];

% reaction time
subplot 121; hold on                                                                % Normal reaction time
plot(1:2,Group_RT(:,1:2),'-o','Color',cc)
plot(1:2,mean(Group_RT(:,1:2)),'*-r','LineWidth',2);                                % auditory congruent (1) incongruent (2)
plot(3:4,Group_RT(:,3:4),'-o','Color',cc)
plot(3:4,mean(Group_RT(:,3:4)),'*-r','LineWidth',2);                                % visual congruent (3) incongruent (4)
xlim([0.5 4.5]); ylim([0.4 0.9])
set(gca,'XTick',1:4,'XTickLabel',{'AUD C','AUD IC','VIS C','VIS IC'});
title('RT with mean');

% performance
subplot 122; hold on;                                                               % Performance
plot(1:2,Group_Perf(:,1:2),'-o','Color',cc)
plot(1:2,median(Group_Perf(:,1:2)),'*-r','LineWidth',2);                            % auditory congruent (1) incongruent (2)
plot(3:4,Group_Perf(:,3:4),'-o','Color',cc)
plot(3:4,median(Group_Perf(:,3:4)),'*-r','LineWidth',2);                            % visual congruent (3) incongruent (4)
xlim([0.5 4.5]); ylim([0.5 1])
set(gca,'XTick',1:4,'XTickLabel',{'AUD C','AUD IC','VIS C','VIS IC'});
title('Perf with Mean');

saveppt('GROUP_CUTOFF_19')


%--------------------------------------------------------------------------
% Difference Plots (Incongruent - Congruent)
%--------------------------------------------------------------------------
Group_DIFF(:,1) = Group_RT(:,2) - Group_RT(:,1);                                    % AIC - AC
Group_DIFF(:,2) = Group_RT(:,4) - Group_RT(:,3);                                    % VIC - VC

figure

% auditory reaction time difference
subplot 211
bar(Group_DIFF(:,1)); xlim([0 21]); set(gca,'XTick',1:20);
title('Mean AIC - Mean AC GROUP')
ylabel('AUDITORY RT DIFFERENCE')
ylim([-0.4 0.4])

% visual reaction time difference
subplot 212
bar(Group_DIFF(:,2)); xlim([0 21]); set(gca,'XTick',1:20);
title('Mean VIC - Mean VC GROUP')
ylabel('VISUAL RT DIFFERENCE')
ylim([-0.4 0.4])

suptitle('NORMAL RT')
saveppt('GROUP_CUTOFF')

%--------------------------------------------------------------------------
% Stats: Reaction times
clc
% signtests on non transformed mean data
[Paud,Haud,STATSaud] = signrank(Group_RT(:,1),Group_RT(:,2));                       % AC vs AIC
Raud = STATSaud.zval / sqrt(40);

[Pvis,Hvis,STATSvis] = signrank(Group_RT(:,3),Group_RT(:,4)) ;                      % VC vs VIC
Rvis = STATSvis.zval / sqrt(40);

%--------------------------------------------------------------------------
% Stats: Performance Scores
clc
% signtests on non transformed median data
[Paud2,Haud2,STATSaud2] = signrank(Group_Perf(:,1),Group_Perf(:,2));                % AC vs AIC
Raud2 = STATSaud2.zval / sqrt(40);

[Pvis2,Hvis2,STATSvis2] = signrank(Group_Perf(:,3),Group_Perf(:,4));                % VC vs VIC
Rvis2 = STATSvis2.zval / sqrt(40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group Decoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear
cd('W:/CIRCLE_EXP/LOG/EEG')
subs = [1:18,20];
n = length(subs);
t = 261; e = 128;                                                                   % t = time, e = electrodes

% % initialise things
SD_AZ = zeros(n,t); SD_A = zeros(e,t,n); SD_Y = cell(1,n);
trials = [];
trialInd = cell(1,n);

for SUBJ = 1:n;
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)));
    load(sprintf('S0%d_DECODE_CUTOFF.mat',subs(SUBJ)))
    
    SD_AZ(SUBJ,:) = SD.Az;                                                          % Aud vs. Vis decoder, Az
    SD_A(:,:,SUBJ) = SD.A;                                                          % Aud vs. Vis decoder, Topo
    SD_Y{SUBJ} = SD.Y;                                                              % Aud vs. Vis decoder, Y signals
    
    
    group_ND.audAz(SUBJ,:) = ND.AudAz;                                              % AH vs. AL decoder, Az
    group_ND.visAZ(SUBJ,:) = ND.VisAz;                                              % VS vs. VL decoder, Az
    group_ND.audTopo(:,:,SUBJ) = ND.TopoAud;                                        % AH vs. AL decoder, Topo
    group_ND.visTopo(:,:,SUBJ) = ND.TopoVis;                                        % VS vs. VL decoder, Topo
    group_ND.audY{SUBJ} = ND.Yaud;                                                  % AH vs. AL decoder, Y signals
    group_ND.visY{SUBJ} = ND.Yvis;                                                  % VS vs. VL decoder, Y signals
    group_ND.YconA{SUBJ} = ND.YConA;
    group_ND.YconV{SUBJ} = ND.YConV;
    group_ND.trials{SUBJ} = trialinfo;
    
    if SUBJ==1;
        group_ND.TAX = Decode_times;
        group_ND.trange = trange;
    end
    
    trials = cat(1,trials,trialinfo);
    trialInd{SUBJ} = trialinfo;
    fprintf('S0%d...\n',SUBJ);
end

cd('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding')

% Group values
G_SD_Y= cat(1,SD_Y{:});                                                             % Group Sensory Decoder (Aud vs. Vis) Y signals
G_YA = cat(1,ND_YA{:});                                                             % Group Normal Decoder (AH vs AL) AUD Y signals
G_YV = cat(1,ND_YV{:});                                                             % Group Normal Decoder (VS vs VL) VIS Y signals

save('GROUP_DECODE_CUTOFF_19.mat','group_ND');

% Auditory trials Y
groupA = []; groupA_Ycon = []; aTrials = [];
groupA = cat(1,groupA,group_ND.audY{:});
groupA_Ycon = cat(1,groupA_Ycon,group_ND.YconA{:});
for s = 1:n;
    aTrials = cat(1,aTrials,group_ND.trials{s}(group_ND.trials{s}(:,2)<=2,:));
end

% Visual trials Y
groupV = []; groupV_Ycon =[];
groupV = cat(1,groupV,group_ND.visY{:});
groupV_Ycon = cat(1,groupV_Ycon,group_ND.YconV{:});
vTrials = [];
for s = 1:n;
    vTrials = cat(1,vTrials,group_ND.trials{s}(group_ND.trials{s}(:,2)>=3,:));
end

group_ND.aTrials = aTrials;                                                         % auditory trials (behaviour)
group_ND.aCon = groupA_Ycon;                                                        % auditory labels (decoder)
group_ND.Ya = groupA;                                                               % auditory signals (decoder)

group_ND.vTrials = vTrials;                                                         % visual trials (behaviour)
group_ND.vCon = groupV_Ycon;                                                        % visual labels (decoder)
group_ND.Yv = groupV;                                                               % visual signals (decoder)

group_ND.cIndA = find(aTrials(:,3)==1);                                             % auditory congruent indices
group_ND.icIndA = find(aTrials(:,3)==2);                                            % auditory incongruent indices
group_ND.cIndV = find(vTrials(:,3)==1);                                             % visual congruent indices
group_ND.icIndV = find(vTrials(:,3)==2);                                            % visual incongruent indices

save('GROUP_DECODE_CUTOFF_19.mat','group_ND','TIMES')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO THE AZ BOOT TO GET THE SIGNIFICANT TIME POINTS
% AZ BOOTSTRAP: compute the Az's 1000 times for data in which the
% condition is randomized across trials. Then compute the group averaged
% for each randomization, take the max-az over time, which yields a
% distribution of 1000 randomized Az values. From these you can compute the
% 1st and 99th percentile and add these to the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ; clc;
addpath('/analyse/Project0109/Lab/ckmatlab/eegck')
addpath('/analyse/Project0109/Lab/ckmatlab/ckinfo')
cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG')
load('GROUP_DECODE_CUTOFF.mat')                                                     % load group Y data

% initialise stuff
btsp = 1000;                                                                        % how many times
n = size(group_ND.audY,2);
T = size(group_ND.TAX,2);

% AZ Randomisations
for subj =1:n;
    
    AzBootA = zeros(btsp,T);
    AzBootV = zeros(btsp,T);
    
    %----------------------------------------------------------------------
    % Auditory bootstrap
    Y = group_ND.audY{subj};                                                        % subject decoding signal from diagonal decoder
    YCon = group_ND.YconA{subj};
    Ntrl = size(YCon,1);
    for b = 1:btsp;                                                                 % each bootstrap
        order = randperm(Ntrl);                                                     % create random order
        YCon2 = YCon(order);                                                        % shuffle labels
        for tw = 1:T                                                                % each decoding timepoint
            [~,~,tmp,~] = ck_decoding_roc(YCon2,Y(:,tw));
            AzBootA(b,tw) = tmp;
        end
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Visual bootstrap
    Y = group_ND.visY{subj};
    YCon = group_ND.YconV{subj};
    Ntrl = size(YCon,1);
    
    for b = 1:btsp;                                                                 % each bootstrap
        order = randperm(Ntrl);                                                     % create random order
        YCon2 = YCon(order);                                                        % shuffle labels
        for tw = 1:T                                                                % each decoding timepoint
            [~,~,tmp,~] = ck_decoding_roc(YCon2,Y(:,tw));
            AzBootV(b,tw) = tmp;
        end
    end
    %----------------------------------------------------------------------
    
    fprintf('Done S0%d Bootstrap Main \n',subj);
    save(sprintf('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_AzBoot/azBoot_S0%d.mat',subj),'AzBootA','AzBootV');
end

%--------------------------------------------------------------------------
% Make a group file
subs = [1:18,20];
n = length(subs);
cd('W:/CIRCLE_EXP/LOG/EEG/Files_AzBoot/')
files = dir('azBoot_S0*.mat');
audAz = zeros(1000,261, n);
visAz = zeros(1000,261, n);

for subj = 1:n;
    load(sprintf('azBoot_S0%d',subs(subj)));
    audAz(:,:,subj) = AzBootA;
    visAz(:,:,subj) = AzBootV;
    fprintf('S0%d done Bootstrap Group \n',subs(subj))
end

%-----------------------------------------------------
% Group Average
azGroupA = sq(mean(audAz,3));
azGroupV = sq(mean(visAz,3));
aMax = max(azGroupA(:,:,1),[],2);
vMax = max(azGroupV(:,:,1),[],2);                                                   % max over time for each randomisation
prcntAUD = prctile(aMax,[1 99]);
prcntVIS = prctile(vMax,[1 99]);

% Plot it
clf
subplot 121; hold on
histfit(aMax); hold on
vline(prcntAUD);
xlabel(sprintf('%s   %s',num2str(prcntAUD(1)),num2str(prcntAUD(2))),'Color','r')
title('AZ BOOT GROUP AUD DISTRIBUTION');

subplot 122; hold on
histfit(vMax)
vline(prcntVIS);
xlabel(sprintf('%s   %s',num2str(prcntVIS(1)),num2str(prcntVIS(2))),'Color','r')
title('AZ BOOT GROUP VIS DISTRIBUTION');

saveas(gcf,'Figure_AzBoot_percentiles.fig');
save('Group_Az_19.mat','audAz','visAz','azGroupA','azGroupV','aMax','vMax','prcntAUD','prcntVIS');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GROUP Y SIGNALS AND SPLIT IT UP BY CONDITION/CONGRUENCY : AUDITORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
cd('W:/CIRCLE_EXP/LOG/EEG');
subs = [1:18, 20];
n = length(subs);
addpath('W:\CIRCLE_EXP\LOG\EEG')
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')
HighToneG = cell(1,n); LowToneG = cell(1,n); CongruentG = cell(1,n); IncongruentG = cell(1,n);
HighToneFG = cell(1,n); LowToneFG = cell(1,n); CongruentFG = cell(1,n); IncongruentFG = cell(1,n);

for SUBJ = 1:n
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    load(sprintf('S0%d_DECODE_CUTOFF.mat',subs(SUBJ)));
    
    % Get y signals from the normal decoder (AH vs AL)
    TIMES = ND.Decode_times;
    audT = find(trialinfo(:,2)<=2);                                                 % auditory trials
    
    % Flip the Y signals that were assigned to a negative value (YconA==0)
    indA = find(ND.YConA==0);
    Y_AUD = ND.Yaud;
    YflipAUD = ND.Yaud;
    YflipAUD(indA,:) = YflipAUD(indA,:)*-1;
    
    % normal trials
    HighTone = Y_AUD(ND.YConA==1,:);                                                % auditory high tone
    LowTone = Y_AUD(ND.YConA==0,:);                                                 % auditory low tone
    audCongruent = Y_AUD(trialinfo(audT,3)==1,:);                                   % auditory congruent
    audIncongruent = Y_AUD(trialinfo(audT,3)==2,:);                                 % auditory incongruent
    
    % flipped trials
    HighToneF = YflipAUD(ND.YConA==0,:);
    LowToneF = YflipAUD(ND.YConA==1,:);
    CongruentF = YflipAUD(trialinfo(audT,3)==1,:);
    IncongruentF = YflipAUD(trialinfo(audT,3)==2,:);
    
    % group files : normal Y signals
    HighToneG{SUBJ} = HighTone;
    LowToneG{SUBJ} = LowTone;
    CongruentG{SUBJ} = audCongruent;
    IncongruentG{SUBJ} = audIncongruent;
    
    % group files : flipped Y signals
    HighToneFG{SUBJ} = HighToneF;
    LowToneFG{SUBJ} = LowToneF;
    CongruentFG{SUBJ} = CongruentF;
    IncongruentFG{SUBJ} = IncongruentF;
    
    suptitle(sprintf('S0%d AUD TRIALS',SUBJ))
    
    cd('W:/CIRCLE_EXP/LOG/EEG/')
    fprintf('S0%d done \n',subs(SUBJ))
end

% get the average of all trials for each subject
meanHT = zeros(20,261); meanLT = meanHT; meanC = meanHT; meanIC = meanHT;
meanHTflip = meanHT; meanLTflip = meanLT; meanCflip = meanHT; meanICflip = meanHT;
for k = 1:n;
    meanHT(k,:) = mean(HighToneG{k});                                               % group mean High Tone
    meanLT(k,:) = mean(LowToneG{k});                                                % group mean Low Tone
    meanC(k,:) = mean(CongruentG{k});                                               % group mean Congruent
    meanIC(k,:) = mean(IncongruentG{k});                                            % group mean Incongruent
    
    meanHTflip(k,:) = mean(HighToneFG{k});                                          % group mean High Tone flipped
    meanLTflip(k,:) = mean(LowToneFG{k});                                           % group mean Low Tone flipped
    meanCflip(k,:) = mean(CongruentFG{k});                                          % group mean Congruent flipped
    meanICflip(k,:) = mean(IncongruentFG{k});                                       % group mean Incongruent flipped
end
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding')

save('Group_Ysignals_19.mat','HighToneG','LowToneG','CongruentG','IncongruentG',...
    'HighToneFG','LowToneFG','CongruentFG','IncongruentFG',...
    'meanHT','meanLT','meanC','meanIC','meanHTflip','meanLTflip','meanCflip','meanICflip')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GROUP Y SIGNALS AND SPLIT IT UP BY CONDITION/CONGRUENCY : VISUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
cd('W:/CIRCLE_EXP/LOG/EEG');
subs = [1:18,20];
n = length(subs);
addpath('W:\CIRCLE_EXP\LOG\EEG')
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')

LargeCircleG = cell(1,n); SmallCircleG = cell(1,n); CongruentG = cell(1,n); IncongruentG = cell(1,n);
LargeCircleFG = cell(1,n); SmallCircleFG = cell(1,n); CongruentFG = cell(1,n); IncongruentFG = cell(1,n);

for SUBJ = 1:n
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    load(sprintf('S0%d_DECODE_CUTOFF.mat',subs(SUBJ)));
    
    % Get y signals from the normal decoder (AH vs AL)
    TIMES = ND.Decode_times;
    visT = find(trialinfo(:,2)>=3);                                                 % visual trials
    
    % Flip the Y signals that were assigned to a negative value
    indV = find(ND.YConV==0);
    Y_VIS = ND.Yvis;
    YflipVIS = ND.Yvis;
    YflipVIS(indV,:) = YflipVIS(indV,:)*-1;
    
    % normal trials
    LargeC = Y_VIS(ND.YConV==0,:);                                                  % visual large circle
    SmallC = Y_VIS(ND.YConV==1,:);                                                  % visual small circle
    visCongruent = Y_VIS(trialinfo(visT,3)==1,:);                                   % visual congruent
    visIncongruent = Y_VIS(trialinfo(visT,3)==2,:);                                 % visual incongruent
    
    % flipped trials
    LargeCF = YflipVIS(ND.YConV==0,:);
    SmallCF = YflipVIS(ND.YConV==1,:);
    CongruentF = YflipVIS(trialinfo(visT,3)==1,:);
    IncongruentF = YflipVIS(trialinfo(visT,3)==2,:);
    
    % group files : normal Y signals
    LargeCircleG{SUBJ} = LargeC;
    SmallCircleG{SUBJ} = SmallC;
    CongruentG{SUBJ} = visCongruent;
    IncongruentG{SUBJ} = visIncongruent;
    
    % group files : flipped Y signals
    LargeCircleFG{SUBJ} = LargeCF;
    SmallCircleFG{SUBJ} = SmallCF;
    CongruentFG{SUBJ} = CongruentF;
    IncongruentFG{SUBJ} = IncongruentF;
    
    cd('W:/CIRCLE_EXP/LOG/EEG/')
    %     saveppt('Y TESTING')
    fprintf('Subj %d done \n',subs(SUBJ));
    
end

%--------------------------------------------------------------------------
% get the average of all trials for each subject
meanLC = zeros(n,261); meanSC = meanLC; meanC = meanLC; meanIC = meanLC;
meanLCflip = meanLC; meanSCflip = meanSC; meanCflip = meanLC; meanICflip = meanLC;
for k = 1:n;
    meanLC(k,:) = mean(LargeCircleG{k});                                            % group mean High Tone
    meanSC(k,:) = mean(SmallCircleG{k});                                            % group mean Low Tone
    meanC(k,:) = mean(CongruentG{k});                                               % group mean Congruent
    meanIC(k,:) = mean(IncongruentG{k});                                            % group mean Incongruent
    
    meanLCflip(k,:) = mean(LargeCircleFG{k});                                       % group mean High Tone flipped
    meanSCflip(k,:) = mean(SmallCircleFG{k});                                       % group mean Low Tone flipped
    meanCflip(k,:) = mean(CongruentFG{k});                                          % group mean Congruent flipped
    meanICflip(k,:) = mean(IncongruentFG{k});                                       % group mean Incongruent flipped
end

save('Group_Y_VIS_19.mat','LargeCircleG','SmallCircleG','CongruentG','IncongruentG',...
    'LargeCircleFG','SmallCircleFG','CongruentFG','IncongruentFG',...
    'meanLC','meanSC','meanC','meanIC','meanLCflip','meanSCflip','meanCflip','meanICflip','TIMES')
%--------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GROUP REGRESSION WEIGHTS, Y signals, and signal differences : AUDITORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG');
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')
subs = [1:18,20];
n = length(subs);

% initialise things
group_AW = zeros(n,261); group_AWD = zeros(n,261,2); group_AWC = zeros(n,261); group_AWIC = zeros(n,261);
group_Yaud = cell(1,n); group_Yflip = cell(1,n);
group_ACind = cell(1,n); group_AICind = cell(1,n);
group_aDiff = zeros(1,n); group_ayDiff = zeros(n,261); group_ayFlip = zeros(n,261);
group_aDiff2 = zeros(1,n); group_ayDiff2 = zeros(n,261);group_ayFlip2 = zeros(n,261);
group_ASTATS = zeros(n,261,4); group_ASTATS_Dual = group_ASTATS;

% Loop through subjects
for SUBJ = 1:n;
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    
    if SUBJ ==1;
        load('S01_DECODE_CUTOFF.mat','Decode_times')                                % load the timings
    end
    
    load(sprintf('S0%d_REG_NEW',subs(SUBJ)))                                        % regression weights
    load(sprintf('S0%d_signalDiffs',subs(SUBJ)))                                    % signal differences
    load(sprintf('S0%d_SignalDiff_btsp',subs(SUBJ)))                                % signal differences : bootstrapped version
    
    % Regression Weights
    group_AW(SUBJ,:) = AUD_Weight(:,1)';                                            % auditory weight, 1 model
    group_AWD(SUBJ,:,1) = AUD_Weight_Dual(:,1)';                                    % auditory weight congruent, dual model
    group_AWD(SUBJ,:,2) = AUD_Weight_Dual(:,2)';                                    % auditory weight incongruent, dual model
    
    % Regression Predictive Power
    group_ASTATS(SUBJ,:,:) = AUD_STATS;                                             % one model
    group_ASTATS_Dual(SUBJ,:,:) = AUD_STATS_Dual;                                   % dual model
    
    % Y signals
    group_Yaud{SUBJ} = Y_AUD;                                                       % Y signals auditory
    group_Yflip{SUBJ} = YflipAUD;                                                   % Y signals auditory Flipped
    
    % Indices
    group_ACind{SUBJ} = AC_ind;                                                     % trial index, congruent auditory
    group_AICind{SUBJ} = AIC_ind;                                                   % trial index, incongruent auditory
    
    % Signal Differences
    group_aDiff(SUBJ) = aDiff;                                                      % reaction time (Aud Incongruent - Aud Congruent)
    group_ayDiff(SUBJ,:) = ayDiff;                                                  % Y signal diff (Aud Incongruent - Aud Congruent)
    group_ayFlip(SUBJ,:) = ayFlipDiff;                                              % Y signal diff Flipped (Aud Incongruent - Aud Congruent)
    
    % Signal Differences Bootstrap
    group_aDiff2(SUBJ) = mean(ADiff);                                               % reaction time (Aud Incongruent - Aud Congruent)
    group_ayDiff2(SUBJ,:) = mean(AY_Diff);                                          % Y signal diff (Aud Incongruent - Aud Congruent)
    group_ayFlip2(SUBJ,:) = mean(AY_DiffF);                                         % Y signal diff Flipped (Aud Incongruent - Aud Congruent)
    
    fprintf('S0%d done \n',subs(SUBJ))
end

cd('W:\CIRCLE_EXP\LOG\EEG');
save('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Regression_19.mat','group_AW','group_AWD','Decode_times','group_ASTATS','group_ASTATS_Dual')
save('W:/CIRCLE_EXP/LOG/EEG/Files_Decoding/Group_SignalDiffs_19.mat','group_aDiff','group_ayDiff','group_ayFlip','group_aDiff2','group_ayDiff2','group_ayFlip2','Decode_times');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GROUP REGRESSION WEIGHTS, Y signals, and signal differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG');
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')
subs = [1:18,20];
n = length(subs);

% initialise things
group_VW = zeros(n,261); group_VWD = zeros(n,261,2);
group_Yvis = cell(1,n); group_Yvisflip = cell(1,n);
group_VCind = cell(1,n); group_VICind = cell(1,n);
group_vDiff = zeros(1,n); group_vyDiff = zeros(n,261); group_vyFlip = zeros(n,261);
group_vDiff2 = zeros(1,n); group_vyDiff2 = zeros(n,261);group_vyFlip2 = zeros(n,261);
group_VSTATS = zeros(n,261,4); group_VSTATS_Dual = group_VSTATS;

% Regression models:
% 1 model | where Y activity auditory is regressed against Reaction time.
% Includes 1 preditor (RT) and constant (1).
% Dual model | where Y activity auditory is regressed against Reaction time
% split up into congruent and incongruent conditions.
% Includes 2 predictors (RT congruent | RT incongruent) and constant (1).

for SUBJ = 1:n
    
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    
    if SUBJ ==1;
        load('S01_DECODE_CUTOFF.mat','Decode_times')                                % load the timings
    end
    
    load(sprintf('S0%d_REG_VIS',subs(SUBJ)))                                        % regression weights
    load(sprintf('S0%d_signalDiffs_VIS',subs(SUBJ)))                                % signal differences
    %     load(sprintf('S0%d_SignalDiff_btsp',SUBJ))                                    % signal differences : bootstrapped version
    
    % Regression Weights
    group_VW(SUBJ,:) = VIS_Weight(:,1)';                                            % auditory weight, 1 model
    group_VWD(SUBJ,:,1) = VIS_Weight_Dual(:,1)';                                    % auditory weight congruent, dual model
    group_VWD(SUBJ,:,2) = VIS_Weight_Dual(:,2)';                                    % auditory weight incongruent, dual model
    
    % Regression Predictive Power
    group_VSTATS(SUBJ,:,:) = VIS_STATS;                                             % one model
    group_VSTATS_Dual(SUBJ,:,:) = VIS_STATS_Dual;                                   % dual model
    
    % Y signals
    group_Yvis{SUBJ} = Y_VIS;                                                       % Y signals auditory
    group_Yvisflip{SUBJ} = YflipVIS;                                                % Y signals auditory Flipped
    
    % Indices
    group_VCind{SUBJ} = VC_ind;                                                     % trial index, congruent auditory
    group_VICind{SUBJ} = VIC_ind;                                                   % trial index, incongruent auditory
    
    % Signal Differences
    group_vDiff(SUBJ) = vDiff;                                                      % reaction time (Aud Incongruent - Aud Congruent)
    group_vyDiff(SUBJ,:) = vyDiff;                                                  % Y signal diff (Aud Incongruent - Aud Congruent)
    group_vyFlip(SUBJ,:) = vyFlipDiff;                                              % Y signal diff Flipped (Aud Incongruent - Aud Congruent)
    
    fprintf('S0%d done \n',subs(SUBJ))
end

cd('W:\CIRCLE_EXP\LOG\EEG');
save('Group_Regression_VIS_19.mat','group_VW','group_VWD','Decode_times','group_VSTATS','group_VSTATS_Dual')
save('Group_SignalDiffs_VIS_19.mat','group_vDiff','group_vyDiff','group_vyFlip','Decode_times');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y FLIPPED Signals: Correct for Multiple Comparisons : AUDITORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding');
load('Group_Ysignals_19.mat','meanCflip','meanICflip')
load('Group_Ysignals_19.mat','TIMES')

I = find(TIMES>=0 & TIMES<=0.545);
TAX = TIMES(I);

n = size(meanCflip,1);
TshufA =zeros(1,length(I),1000); TtrueA=zeros(1,length(I));

audCongruent = meanCflip(:,I);                                                      % group Congruent   Y signals
audIncongruent = meanICflip(:,I);                                                   % group Incongruent Y signals

for k=1:1000
    
    for t=1:length(I)
        
        % Work out the true tvalue
        if k==1
            [~,~,~,STATS] = ttest(audCongruent(:,t),audIncongruent(:,t));
            TtrueA(1,t) = STATS.tstat;
        end
        
        % avoid the tru order of subjects
        order = randperm(n);
        if sum(abs(order-(1:n)))==0
            order = randperm(n);
        end
        
        % Work out the shuffled tvalue (shuffle Y by subject)
        [~,~,~,STATS2] = ttest(audCongruent(order,t),audIncongruent(:,t));
        TshufA(1,t,k) = STATS2.tstat; % has to be a 3D matrix
    end
end

cfg.critval = 1.8;                                                                  % threshold for signif t-values
cfg.clusterstatistic = 'maxsum';                                                    % maxsize maxsum
cfg.critvaltype = 'par';                                                            % parametric threshold
cfg.minsize = 2;                                                                    % minimum cluster size
cfg.pval = 0.05;                                                                    % threshold to select signifciant clusters
cfg.df = 19;
[audPC_Y,audNC_Y] = eegck_clusterstats(cfg,TtrueA,TshufA);

audCongruentM = mean(audCongruent);
audIncongruentM = mean(audIncongruent);

% Plot it
clf; hold on;
suptitle('Group Y auditory signals, corrected for multiple comparisons')
plot(TAX,audCongruentM,'b','LineWidth',2);
plot(TAX,audIncongruentM,'r','LineWidth',2);
ind1 = find(audPC_Y.maskSig~=0);
plot(TAX(ind1),audCongruentM(ind1),'k*');
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
hline(0); vline(0); vline(0.3);
% saveppt('group_new')

audPC_Y.ind = ind1; 

save('Group_Shuffled_Ysignals_19_new.mat','TshufA','TtrueA','audPC_Y','audNC_Y','cfg',...
    'audCongruent','audIncongruent','audCongruentM','audIncongruentM','I','TAX')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y FLIPPED Signals: Correct for Multiple Comparisons : VISUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding');
load('Group_Ysignals_VIS_19.mat','meanCflip','meanICflip')
load('TIMES');

I = find(TIMES>=0 & TIMES<=0.545);
TAX = TIMES(I);
n = size(meanCflip,1);
TshufV =zeros(1,length(I),1000); TtrueV=zeros(1,length(I));

visCongruent = meanCflip(:,I);                                                      % group Congruent   Y signals
visIncongruent = meanICflip(:,I);                                                   % group Incongruent Y signals

for k=1:1000
    
    for t=1:length(I)
        
        % Work out the true tvalue
        if k==1
            [~,~,~,STATS] = ttest(visCongruent(:,t),visIncongruent(:,t));
            TtrueV(1,t) = STATS.tstat;
        end
        
        % avoid the tru order of subjects
        order = randperm(n);
        if sum(abs(order-(1:n)))==0
            order = randperm(n);
        end
        
        % Work out the shuffled tvalue (shuffle Y by subject)
        [~,~,~,STATS2] = ttest(visCongruent(order,t),visIncongruent(:,t));
        TshufV(1,t,k) = STATS2.tstat;                                               % has to be a 3D matrix
    end
end

cfg.critval = 1.8;                                                                  % threshold for signif t-values
cfg.clusterstatistic = 'maxsum';                                                    % maxsize maxsum
cfg.critvaltype = 'par';                                                            % parametric threshold
cfg.minsize = 2;                                                                    % minimum cluster size
cfg.pval = 0.05;                                                                    % threshold to select signifciant clusters
cfg.df = 19;
[visPC_Y,visNC_Y] = eegck_clusterstats(cfg,TtrueV,TshufV)

visCongruentM = mean(visCongruent);
visIncongruentM = mean(visIncongruent);

% Plot it
clf; hold on;
suptitle('Group Y visual signals, corrected for multiple comparisons')
plot(TAX,visCongruentM,'b','LineWidth',2);
plot(TAX,visIncongruentM,'r','LineWidth',2);
ind1 = find(visPC_Y.maskSig~=0);
ind2 = find(visNC_Y.maskSig~=0);
plot(TAX(ind1),visCongruentM(ind1),'k*');
plot(TAX(ind2),visCongruentM(ind2),'k*');
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
hline(0); vline(0); vline(0.3);
% saveppt('group_new')

visPC_Y.ind = ind1; 
visNC_Y.ind = ind2; 


save('Group_Shuffled_Ysignals_VIS_19_new.mat','TshufV','TtrueV','visPC_Y','visNC_Y',...
    'cfg','visCongruent','visIncongruent','visCongruentM','visIncongruentM','I','TAX')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression: Correct for Multiple Comparisons : AUDITORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Regression');
subs = [1:18,20];
n = length(subs);
Group_Aboot = zeros(n,261,1000);
Group_ACboot = Group_Aboot;
Group_AICboot = Group_Aboot;

% Get the group shuffled regression weight matrices
for SUBJ = 1:n;
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    load(sprintf('S0%d_Regression_Boot',subs(SUBJ)))
    
    Group_Aboot(SUBJ,:,:)   = Aboot;                                            % one model, auditory shuffled weight
%     Group_ACboot(SUBJ,:,:)  = ACboot;                                           % dual model, auditory congruent shuffled weight
%     Group_AICboot(SUBJ,:,:) = AICboot;                                          % dual model, auditory incongruent shuffled weight
%     
    fprintf('S0%d done \n',subs(SUBJ))
end

% Generate the tvalues from the random regression weights
cd('W:/CIRCLE_EXP/LOG/EEG');
Tshuf = zeros(1000,length(I)); Tshuf_C = Tshuf; TshufC = Tshuf; TshufIC = Tshuf;
for b = 1:1000;
    for t = 1:length(I);
        [~,~,~,stats] = ttest(Group_Aboot(:,I(t),b));                               % test against chance
        Tshuf(b,t) = stats.tstat;
        
        [~,~,~,stats2] = ttest(Group_ACboot(:,I(t),b),Group_AICboot(:,I(t),b));     % test against each other
        Tshuf_C(b,t) = stats2.tstat;
    end
end
% ~10mins

save('Group_Shuffled_RegWeights_19.mat',...
    'Group_Aboot','Group_ACboot','Group_AICboot',...
    'Tshuf','Tshuf_C',...
    'TAX','I')


% load in the real regression weights and get the real tvalues
load('Group_Regression_19.mat','group_AW','group_AWD')
Ttrue = zeros(1,length(I)); TtrueC = zeros(1,length(I)); TtrueIC= Ttrue; Ttrue_C = Ttrue;
for t = 1:length(I);
    
    % neural regression weights
    [~,~,~,statsT] = ttest(group_AW(:,I(t)));                                       % test against chance
    Ttrue(t) = statsT.tstat;
    
%     [~,~,~,statsT2] = ttest(group_AWD(:,I(t),1),group_AWD(:,I(t),2));               % test against each other
%     Ttrue_C(t) = statsT2.tstat;
    
end

save('Group_Shuffled_RegWeights_19.mat',...
    'Group_Aboot','Group_ACboot','Group_AICboot',...
    'Tshuf','Tshuf_C','TshufC','TshufIC','TAX','I',...
    'Ttrue','TtrueC','TtrueIC','Ttrue_C')


% Test the randomly shuffed distribution against the real one
cfg.critval = 1.8;                                                                  % threshold for signif t-values
cfg.clusterstatistic = 'maxsum';                                                    % maxsize maxsum
cfg.critvaltype = 'par';                                                            % parametric threshold
cfg.minsize = 2;                                                                    % minimum cluster size
cfg.pval = 0.05;                                                                    % threshold to select signifciant clusters
cfg.df = 19;

% one model weights tested against chance
[PC_R,NC_R] = eegck_clusterstats(cfg,Ttrue,Tshuf);                                  % one model weights
ind_OneModel = find(PC_R.maskSig~=0);

% save : one model results
oneModel = [];
oneModel.PC = PC_R;
oneModel.NC = NC_R;
oneModel.Ttrue = Ttrue;
oneModel.Tshuf = Tshuf;
oneModel.ind = ind_OneModel;


save('Group_Shuffled_RegWeights_19_new.mat','Group_Aboot','Group_ACboot','Group_AICboot','TAX','I',...
    'oneModel','Ttrue','Tshuf')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression: Correct for Multiple Comparisons : VISUAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Regression');
subs = [1:18,20];
n = length(subs);
addpath('W:\CIRCLE_EXP\LOG\EEG')
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')

Group_Vboot = zeros(20,261,1000);
Group_VCboot = Group_Vboot;
Group_VICboot = Group_Vboot;

% Get the group shuffled regression weight matrices
for SUBJ = 1:n;
    cd(sprintf('W:/CIRCLE_EXP/LOG/EEG/Subj_%d',subs(SUBJ)))
    load(sprintf('S0%d_Regression_Boot_VIS',subs(SUBJ)))
    
    Group_Vboot(SUBJ,:,:)   = Vboot;                                            % one model, auditory shuffled weight
    Group_VCboot(SUBJ,:,:)  = VCboot;                                           % dual model, auditory congruent shuffled weight
    Group_VICboot(SUBJ,:,:) = VICboot;                                          % dual model, auditory incongruent shuffled weight
    
    fprintf('S0%d done \n',subs(SUBJ))
end


% Generate the tvalues from the random regression weights
cd('W:/CIRCLE_EXP/LOG/EEG');
load('Group_Correlations_VIS_19.mat','TAX');

I = find(TAX>=0 & TAX<=0.545);
TshufV = zeros(1000,length(I)); Tshuf_VC = TshufV; TshufVC = TshufV; TshufVIC = TshufV;
for b = 1:1000;
    for t = 1:length(I);
        [~,~,~,stats] = ttest(Group_Vboot(:,I(t),b));                           % test against chance
        TshufV(b,t) = stats.tstat;
        
        [~,~,~,stats2] = ttest(Group_VCboot(:,I(t),b),Group_VICboot(:,I(t),b));  % test against each other
        Tshuf_VC(b,t) = stats2.tstat;
    end
end
% ~10mins



% load in the real regression weights and get the real tvalues
load('Group_Regression_VIS_19.mat','group_VW','group_VWD')
TtrueV = zeros(1,length(I)); TtrueVC = zeros(1,length(I)); TtrueVIC= TtrueV; Ttrue_VC = TtrueV;
for t = 1:length(I);
    
    [~,~,~,statsT] = ttest(group_VW(:,t));                                          % test against chance
    TtrueV(t) = statsT.tstat;
    
    [~,~,~,statsT2] = ttest(group_VWD(:,t,1),group_VWD(:,t,2));                     % test against each other
    Ttrue_VC(t) = statsT2.tstat;
    
end

save('Group_Shuffled_RegWeights_VIS_19.mat',...
    'Group_Vboot','Group_VCboot','Group_VICboot',...
    'TshufV','Tshuf_VC','TAX','I',...
    'TtrueV','Ttrue_VC')


% Test the randomly shuffed distribution against the real one
cfg.critval = 1.8;                                                                  % threshold for signif t-values
cfg.clusterstatistic = 'maxsum';                                                    % maxsize maxsum
cfg.critvaltype = 'par';                                                            % parametric threshold
cfg.minsize = 2;                                                                    % minimum cluster size
cfg.pval = 0.05;                                                                    % threshold to select signifciant clusters
cfg.df = 19;

% one model weights tested against chance
[PC_R,NC_R] = eegck_clusterstats(cfg,TtrueV,TshufV);                                % one model weights
ind_OneModel = find(PC_R.maskSig~=0 | NC_R.maskSig~=0);

% dual model weights tested against each other (paired ttest)
[PC_R_Dual,NC_R_Dual] = eegck_clusterstats(cfg,Ttrue_VC,Tshuf_VC);
ind_DualModel = find(PC_R_Dual.maskSig~=0 | NC_R_Dual.maskSig~=0);

% save : one model results
oneModel.PC = PC_R;
oneModel.NC = NC_R;
oneModel.Ttrue = TtrueV;
oneModel.Tshuf = TshufV;
oneModel.ind = 'NA';

% save: dual model results
dualModel.PC = PC_R_Dual;
dualModel.NC = NC_R_Dual;
dualModel.Ttrue = Ttrue_VC;
dualModel.Tshuf = Tshuf_VC;
dualModel.ind = ind_DualModel;


save('Group_Shuffled_RegWeights_VIS_19_new.mat','Group_Vboot','Group_VCboot','Group_VICboot','TAX','I',...
    'oneModel','dualModel','TtrueV','Ttrue_VC')





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MUTUAL INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% needs done on the grid. 

clear all; close all; clc;
addpath('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/gcmi-master/matlab')
cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient')

subs = [1:18,20];
n = length(subs);
groupmiAC = zeros(128,110,19);
groupmiAIC = zeros(128,110,19);


for SUBJ = 1:n;
    
    %----------------------------------------------------------------------
    % Do the auditory trials first : calculate MI
    %----------------------------------------------------------------------
    fname = sprintf('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_DataAvgd/S0%d_EEG_matrix_AUD.mat',subs(SUBJ));
    % fname = sprintf('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_DataAvgd/S0%d_EEG_matrix_AUD.mat',subs(SUBJ));
    load(fname);
    I = find(TAX>=0 & TAX<=0.545);
    Nsens = 128; Ntime = length(I);
    
    x = dataAud(:,:,I);                                     % eeg data
    x = permute(x,[2 3 1]);                                 % [trials,time,sensors];
    
    % split the data into different congruency conditions
    dAC = x(aTrials(:,3)==1,:,:);
    dAIC = x(aTrials(:,3)==2,:,:);
    
    gAC = gradient_dim1(permute(dAC,[2 1 3]));
    gAC = permute(gAC,[2 1 3]);
    
    gAIC = gradient_dim1(permute(dAIC,[2 1 3]));
    gAIC = permute(gAIC,[2 1 3]);
    
    % copnorm eeg voltage and gradient separately
    dataAC = copnorm(dAC);
    dataAIC = copnorm(dAIC);
    gradAC = copnorm(gAC);
    gradAIC = copnorm(gAIC);
    
    
    % stimuli labels
    stimAC = aTrials(aTrials(:,3)==1,2)==1;                                     % aud congruent trials, stimulus class lables: 1 = high tone, 0 = low tone
    stimAIC = aTrials(aTrials(:,3)==2,2)==1;
    
    
    % now wherever you do info calculation replace cx with [cx cdx] (concatenating to make a 2d response)
    % so something like
    miAC = zeros(Nsens,Ntime);
    miAIC = miAC;
    for si=1:Nsens
        for ti=1:Ntime
            miAC(si,ti) = mi_model_gd( [dataAC(:,ti,si) gradAC(:,ti,si)], stimAC, 2 );
            miAIC(si,ti) = mi_model_gd( [dataAIC(:,ti,si) gradAIC(:,ti,si)], stimAIC, 2 );
        end
    end
    
    
    save(sprintf('S0%d_MI_gradient.mat',SUBJ),'miAC','miAIC','I','TAX')
    
    
    
    % groupMI(:,:,SUBJ) = MI;
    fprintf('SO%d Done \n',SUBJ);
    
    groupmiAC(:,:,SUBJ) = miAC;
    groupmiAIC(:,:,SUBJ) = miAIC;
    
end

save('groupMI.mat','groupmiAC','groupmiAIC','TAX','I')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIGNIFICANCE TESTING FOR ALL OF THE ABOVE COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
addpath('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/gcmi-master/matlab')
cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient')

clear all; close all; clc;
subs = [1:18,20];
n = length(subs);
groupmiAC = zeros(128,110,19);
groupmiAIC = zeros(128,110,19);

Nperm = 1000;

for SUBJ = 1:length(subs);
    
    %----------------------------------------------------------------------
    % Do the auditory trials first : calculate MI
    %----------------------------------------------------------------------
    fname = sprintf('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_DataAvgd/S0%d_EEG_matrix_AUD.mat',subs(SUBJ));
    load(fname);
    I = find(TAX>=0 & TAX<=0.545);
    Nsens = 128; Ntime = length(I);
    
    x = dataAud(:,:,I);                                     % eeg data
    x = permute(x,[2 3 1]);                                 % [trials,time,sensors];
    
    % split the data into different congruency conditions
    dAC = x(aTrials(:,3)==1,:,:);
    dAIC = x(aTrials(:,3)==2,:,:);
    
    gAC = gradient_dim1(permute(dAC,[2 1 3]));
    gAC = permute(gAC,[2 1 3]);
    
    gAIC = gradient_dim1(permute(dAIC,[2 1 3]));
    gAIC = permute(gAIC,[2 1 3]);
    
    % copnorm eeg voltage and gradient separately
    dataAC = copnorm(dAC);
    dataAIC = copnorm(dAIC);
    gradAC = copnorm(gAC);
    gradAIC = copnorm(gAIC);
    
    
    % stimuli labels
    stimAC = aTrials(aTrials(:,3)==1,2)==1;                                     % aud congruent trials, stimulus class lables: 1 = high tone, 0 = low tone
    stimAIC = aTrials(aTrials(:,3)==2,2)==1;
    
    idxAC = zeros(size(stimAC,1),Nperm);
    idxAIC = zeros(size(stimAIC,1),Nperm);
    pstimAC = idxAC; pstimAIC = idxAIC;
    
    for b = 1:Nperm;
        idxAC(:,b) = randperm(size(stimAC,1));
        idxAIC(:,b) = randperm(size(stimAIC,1));
        pstimAC(:,b) = stimAC(idxAC(:,b));
        pstimAIC(:,b) = stimAIC(idxAIC(:,b));
    end
    
    
    % now wherever you do info calculation replace cx with [cx cdx] (concatenating to make a 2d response)
    % so something like
    permAC = zeros(Nsens,Ntime,Nperm);
    permAIC = permAC;
    
    for b = 1:Nperm;
        for si=1:Nsens
            for ti=1:Ntime
                permAC(si,ti,b) = mi_model_gd( [dataAC(:,ti,si) gradAC(:,ti,si)], pstimAC(:,b), 2 );
                permAIC(si,ti,b) = mi_model_gd( [dataAIC(:,ti,si) gradAIC(:,ti,si)], pstimAIC(:,b), 2 );
            end
        end
    end
    % ~1.5 hours each subject
    
    save(sprintf('S0%d_MI_gradient_perm.mat',SUBJ),'permAC','permAIC','I','TAX')
    
    fprintf('SO%d Done \n',SUBJ);

    
end

% save('groupMI.mat','groupmiAC','groupmiAIC','TAX','I')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE A GROUP FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
cd('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient\Files_IndividualFiles')

groupPermAC = zeros(128,110,1000,19);
groupPermAIC = groupAC;

for SUBJ = 1:19;
    load(sprintf('S0%d_MI_gradient_perm.mat',SUBJ))
    
    groupPermAC(:,:,:,SUBJ) = permAC;
    groupPermAIC(:,:,:,SUBJ) = permAIC;
    
    fprintf('S0%d done \n',SUBJ);
end

cd ..
save('groupPerm.mat','groupPermAC','groupPermAIC','TAX','I','-v7.3')

groupAVG(:,:,1) = mean(groupmiAC,3);            % MI AC
groupAVG(:,:,2) = mean(groupmiAIC,3);           % MI AIC
groupPermAVG(:,:,:,1) = mean(groupPermAC,4);    % MI PERM AC
groupPermAVG(:,:,:,2) = mean(groupPermAIC,4);   % MI PERM AIC

ImaxAC = max(squeeze(max(groupPermAVG(:,:,:,1),[],2)),[],1);
threshA = prctile(ImaxAC, 99.9);                       % 99th percentile
IsigAC = groupAVG(:,:,1)>threshA;

ImaxAIC = max(squeeze(max(groupPermAVG(:,:,:,2),[],2)),[],1);
threshAIC = prctile(ImaxAIC, 99.9);                       % 99th percentile
IsigAIC = groupAVG(:,:,2)>threshAIC;

clf
subplot 221; imagesc(TAX(I),[],groupAVG(:,:,1)); title('Auditory Congruent')
subplot 222; imagesc(TAX(I),[],IsigAC); title('Auditory Congruent')
subplot 223; imagesc(TAX(I),[],groupAVG(:,:,2)); title('Auditory Incongruent')
subplot 224; imagesc(TAX(I),[],IsigAIC); title('Auditory Incongruent')

saveas(gcf,'Figure_MI_gradient_Sig.fig')

[electsC,timeC] = find(IsigAC~=0);
electsC = unique(electsC);
timeC = unique(timeC);

[electsIC,timeIC] = find(IsigAIC~=0);
electsIC = unique(electsIC);
timeIC = unique(timeIC);

save('groupMI_gradient_Sig.mat','groupAVG',...
    'IsigAC','IsigAIC','threshA','threshAIC','TAX','I',...
    'electsC','timeC','electsIC','timeIC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE CONGRUENT VS INCONGRUENT, CORRECTED FOR MULTIPLE COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient');
cd('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient');
load groupMI
groupMI(:,:,:,1) = groupmiAC;
groupMI(:,:,:,2) = groupmiAIC;
TtrueA = zeros(128,110);

for e = 1:128;
    for t = 1:length(I);
        [~,~,~,STATS1] = ttest(groupMI(e,t,:,1),groupMI(e,t,:,2));
        TtrueA(e,t) = STATS1.tstat;
        
    end
end

boot = 1000;
TshufA = zeros(128,110,boot);
for b = 1:boot;
    
    for e = 1:128;
        
        for t = 1:length(I);
            
            order = randperm(19);
            
            % Build the shuffled ones
            [~,~,~,STATS1] = ttest(groupMI(e,t,order,1),groupMI(e,t,:,2));
            TshufA(e,t,b) = STATS1.tstat;
            
            
        end
    end
end

save('tvalues.mat','TshufA','TtrueA','groupMI','TAX','I');

cfg.clusterstatistic = 'maxsum';                                                    % maxsize maxsum
cfg.critvaltype = 'par';
cfg.critval = 1.5;
cfg.minsize = 2;                                                                    % minimum cluster size
cfg.pval = 0.05;                                                                    % threshold to select signifciant clusters
cfg.df = 19;
cfg.conn = 6;

[audPC,audNC] = eegck_clusterstats(cfg,TtrueA,TshufA);
aMask = audPC.maskSig+audNC.maskSig;


load('W:/CIRCLE_EXP/LOG/EEG/128labels.mat')

TAX2 = TAX(I);

clear aPos
aPos.ind = find(audPC.p<0.05);
aPos.p = audPC.p(aPos.ind);
aPos.t = audPC.stat(aPos.ind);
aPos.effect = audPC.Effect(aPos.ind);
aPos.maskSig = audPC.maskSig;
aPos.stat = audPC.stat;
for k = 1:length(aPos.ind);
    [test,test2] = find(aPos.maskSig==aPos.ind(k));
    aPos.clusterTime{k} = TAX2(unique(test2));
    aPos.clusterElecs{k} = labels(unique(test));
    aPos.clusterTimeInd{k} = unique(test2);
    aPos.clusterElecInd{k} = unique(test);
end

clear aNeg
aNeg.ind = find(audNC.p<0.05);
aNeg.p = audNC.p(aNeg.ind);
aNeg.t = audNC.stat(aNeg.ind);
aNeg.effect = audNC.Effect(aNeg.ind);
aNeg.maskSig = audNC.maskSig;
aNeg.stat = audNC.stat;
for k = 1:length(aNeg.ind);
    [test,test2] = find(aNeg.maskSig==aNeg.ind(k));
    aNeg.clusterTime{k} = TAX2(unique(test2));
    aNeg.clusterElecs{k} = labels(unique(test));
    aNeg.clusterTimeInd{k} = unique(test2);
    aNeg.clusterElecInd{k} = unique(test);
end

save('MI_gradient_stats_new.mat','groupMI','TAX','I','audNC','audPC','aPos','aNeg','aMask','cfg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groupDiff = sq(mean(groupMI(:,:,:,2)-groupMI(:,:,:,1),3));

clf
subplot 131; imagesc(TAX(I),[],groupDiff,[-0.02 0.02]); colorbar; title('MI Diff')
subplot 132; imagesc(TAX(I),[],aMask); colorbar; title('Sig Clusters')
subplot 133; imagesc(TAX(I),[],groupDiff.*aMask,[-0.02 0.02]); colorbar


figure
% do the topoplots
load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat');

EVP = EvpDummy;
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.colorbar = 'yes';
cfg.zlim = [-0.02 0.02];
cfg.highlight = 'on';
cfg.highlightsymbol = '*';
cfg.highlightcolor = [1 1 1];
clf
for k = 1:length(aPos.ind);
    subplot(2,2,k);
    EVP.avg = mean(groupDiff(:,aPos.clusterTimeInd{k}),2);
    EVP.time = aPos.clusterTime{k}(1);
    cfg.highlightchannel = aPos.clusterElecs{k};
    ft_topoplotER(cfg,EVP)
    title(sprintf('Auditory Positive Cluster %d, %s to %s',k,num2str(aPos.clusterTime{k}(1)),num2str(aPos.clusterTime{k}(end))));
    
end

cfg.zlim = [-0.04 0.04];

subplot 224;
cfg.zlim = [-0.02 0.02];
EVP.avg = mean(groupDiff(:,aNeg.clusterTimeInd{1}),2);
EVP.time = aNeg.clusterTime{1}(1);
EVP.time = aNeg.clusterTime{1}(1);
cfg.highlightchannel = aNeg.clusterElecs{1};
ft_topoplotER(cfg,EVP)
title(sprintf('Negative Cluster %d, %s to %s',k,num2str(aNeg.clusterTime{1}(1)),num2str(aNeg.clusterTime{1}(end))));

suptitle('Auditory Congruent MI - Auditory Incongruent MI')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REGRESSION ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the significant clusters, and extract the EEG data underlying this,
% and then regress the activity at the three clusters against reaction time

cd('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient')
load MI_gradient_stats

n = 19;

e1 = aPos.clusterElecInd{1};                % cluster 1 significant electrodes
e2 = aPos.clusterElecInd{2};                % cluster 2 significant electrodes
e3 = aNeg.clusterElecInd{1};                % cluster 3 significant electrodes 

t1 = aPos.clusterTimeInd{1}(1:10);          % cluster 1 significant time
t2 = aPos.clusterTimeInd{2}(1:10);          % cluster 2 significant time
t3 = aNeg.clusterTimeInd{1};                % cluster 3 significant time 


for SUBJ = 1:n;
    
    clearvars -except aPos aNeg SUBJ n groupB e1 e2 e3 t1 t2 t3
    load(sprintf('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_DataPrepro/S0%d_EEG_matrix.mat',SUBJ))
    
    RT = trials(:,7);
    I = find(TAX>=0 & TAX<=0.545);
    dataSmall = data(:,:,I);
    
    cluster1 = sq(mean(mean(dataSmall(:,e1,t1),2),3));
    cluster2 = sq(mean(mean(dataSmall(:,e2,t2),2),3));
    cluster3 = sq(mean(mean(dataSmall(:,e3,t3),2),3));
    
    clusters = [cluster1 cluster2 cluster3];
    clusters(:,4) = 1;
    
    [BETA,INTERVALS,RESIDUALS,~,STATS] = regress(RT,clusters);
    
    save(sprintf('S0%d_Regress.mat',SUBJ),'BETA','clusters','INTERVALS','RESIDUALS','STATS')
    fprintf('S0%d done \n',SUBJ);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE A GROUP FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groupB = zeros(10,4); 
for s = 1:19;
    load(sprintf('S0%d_Regress',s),'BETA');
    groupB(s,:) = BETA;
end
clearvars -except groupEEG groupRT groupB


for k = 1:3;
    for s = 1:19;        
        plot(k,groupB(s,k),'o'); hold on;
    end
end
xlim([0.5 3.5]);
hold on
boxplot(groupB(:,1:3),'Colors','k'); ylabel('beta weights');

% calculate t values
t(1) = mean(groupB(:,1))./sem(groupB(:,1));             % cluster 1
t(2) = mean(groupB(:,2))./sem(groupB(:,2));             % cluster 2
t(3) = mean(groupB(:,3))./sem(groupB(:,3));             % cluster 3 

% Do the stats 
[H1,P1,CI1,STATS1] = ttest(groupB(:,1));
d1 =  mean(groupB(:,1))/std(groupB(:,1));
r1 = (t(1)^2) / (t(1)^2 + 18);

[H2,P2,CI2,STATS2] = ttest(groupB(:,2));
d2 =  mean(groupB(:,2))/std(groupB(:,2));
r2 = (t(2)^2) / (t(2)^2 + 18);

[H3,P3,CI3,STATS3] = ttest(groupB(:,3));
d3 =  mean(groupB(:,1))/std(groupB(:,3));
r3 = (t(3)^2) / (t(3)^2 + 18);

load('W:/CIRCLE_EXP/LOG/EEG/128labels.mat')

cluster1.t = t(1);
cluster1.p = P1;
cluster1.CI = CI1;
cluster1.sd = STATS1.sd;
cluster1.d = d1;
cluster1.r = r1;
cluster1.timeInd = t1;
cluster1.time = TAX2(t1);
cluster1.elecInd = e1;
cluster1.elects = labels(e1);


cluster2.t = t(2);
cluster2.p = P2;
cluster2.CI = CI2;
cluster2.sd = STATS2.sd;
cluster2.d = d2;
cluster2.r = r2;
cluster2.timeInd = t2;
cluster2.time = TAX2(t2);
cluster2.elecInd = e2;
cluster2.elects = labels(e2);


cluster3.t = t(3);
cluster3.p = P3;
cluster3.CI = CI3;
cluster3.sd = STATS3.sd;
cluster3.d = d3;
cluster3.r = r3;
cluster3.timeInd = t3;
cluster3.time = TAX2(t3);
cluster3.elecInd = e3;
cluster3.elects = labels(e3);


save('group_Regression_New.mat','cluster1','cluster2','cluster3','groupB','TAX2')

%%

% 
% %% GROUP REGRESSION : REGRESS RT ON GROUP DATA
% y2 = groupEEG;  % eeg data
% x = groupRT;        % response
% 
% nobs = size(y2,1);
% nregions = 3;
% 
% % common slope
% X = cell(nobs,1);
% for j=1:nobs
%     X{j} = [eye(nregions), repmat(groupRT(j),nregions,1)];
% end
% 
% [BETA,sig,resid,vars,loglik2]  = mvregress(X,y2);
% cc = [0.5 0.8 1; 1 0.6 0.8; 0.8 1 0.8];
% 
% B2 = [BETA(1:nregions)';repmat(BETA(end),1,nregions)];
% xx = linspace(0.3,1.3)';
% clf
% for k = 1:3;
%     plot(groupRT,groupEEG(:,k),'x','color',cc(k,:)); hold on
% end
% 
% h = plot(xx, [ones(size(xx)),xx]*B2,'-','LineWidth',2);
% legend(h,{sprintf('C1 beta: %s',num2str(BETA(1))),sprintf('C2 beta: %s',num2str(BETA(2))),sprintf('C3 beta: %s',num2str(BETA(3)))},'FontSize',8)
% xlim([0 1.5]); ylim([-20 20])
% title('EEG data at 3 sig clusters regressed against RT')
% xlabel('RT'); ylabel('EEG')






