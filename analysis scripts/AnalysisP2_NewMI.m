%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW MI SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
addpath('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/gcmi-master/matlab')
cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient')

% addpath('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/gcmi-master/matlab')
% cd('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation')

subs = [1:18,20];
n = length(subs);
groupmiAC = zeros(128,110,19);
groupmiAIC = zeros(128,110,19);

%%
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




%% ------------------------------------------------------------------------
% SIGNIFICANCE TESTING FOR ALL OF THE ABOVE COMPARISONS
%--------------------------------------------------------------------------
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
    
    
    
    % groupMI(:,:,SUBJ) = MI;
    fprintf('SO%d Done \n',SUBJ);
%     
%     groupmiAC(:,:,SUBJ) = miAC; 
%     groupmiAIC(:,:,SUBJ) = miAIC;
    
end

% save('groupMI.mat','groupmiAC','groupmiAIC','TAX','I')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% MAKE A GROUP FILE
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
    'IsigAC','IsigAC','threshA','threshAIC','TAX','I')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE CONGRUENT VS INCONGRUENT, CORRECTED FOR MULTIPLE COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('/analyse/Project0146/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Gradient');
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

save('MI_gradient_stats.mat','groupMI','TAX','I','audNC','audPC','aPos','aNeg','aMask','cfg'); 
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
% EVP.avg = mean(groupMI(:,aNeg.clusterTimeInd{1}),2);
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

e1 = aPos.clusterElecInd{1};
e2 = aPos.clusterElecInd{2};
e3 = aNeg.clusterElecInd{1};

t1 = aPos.clusterTimeInd{1}(1:10);
t2 = aPos.clusterTimeInd{2}(1:10); 
t3 = aNeg.clusterTimeInd{1}; 


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
    
%     
%     y2 = clusters;  % eeg data
%     x = RT;        % response
%     
%     nobs = size(data,1);
%     nregions = 3;
%     
%     % common slope
%     X = cell(nobs,1);
%     for j=1:nobs
%         X{j} = [eye(nregions), repmat(RT(j),nregions,1)];
%     end
%     
%     [BETA,sig,resid,vars,loglik2]  = mvregress(X,y2);
%     
    save(sprintf('S0%d_Regress.mat',SUBJ),'BETA','clusters','INTERVALS','RESIDUALS','STATS')
    
    fprintf('S0%d done \n',SUBJ);
end


%% MAKE A GROUP FILE 
groupEEG = []; groupRT = groupEEG;
for s = 1:19; 
    
    load(sprintf('S0%d_Regress',s),'BETA');
%     groupEEG = cat(1,groupEEG,y2);
%     groupRT = cat(1,groupRT,x); 
    groupB(s,:) = BETA; 
end
clearvars -except groupEEG groupRT groupB

% save('groupRegress.mat','groupEEG','groupRT','groupB')


clf
for k = 1:3;
    for s = 1:19;
        
        plot(k,groupB(s,k),'o'); hold on;
    end
end
xlim([0.5 3.5]); 
hold on
boxplot(groupB(:,1:3),'Colors','k'); ylabel('beta weights'); 

% calculate t values
t(1) = mean(groupB(:,1))./sem(groupB(:,1)); 
t(2) = mean(groupB(:,2))./sem(groupB(:,2)); 
t(3) = mean(groupB(:,3))./sem(groupB(:,3));



[H1,P1,CI1,STATS1] = ttest(groupB(:,1));
d1 =  mean(groupB(:,1))/std(groupB(:,1)); 
r1 = (t(1)^2) / (t(1)^2 + 18);

[H2,P2,CI2,STATS2] = ttest(groupB(:,2));
d2 =  mean(groupB(:,2))/std(groupB(:,2)); 
r2 = (t(2)^2) / (t(2)^2 + 18);

[H3,P3,CI3,STATS3] = ttest(groupB(:,3));
d3 =  mean(groupB(:,3))/std(groupB(:,3)); 
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


%% GROUP REGRESSION : REGRESS RT ON GROUP DATA 
y2 = groupEEG;  % eeg data
x = groupRT;        % response

nobs = size(y2,1);
nregions = 3;

% common slope
X = cell(nobs,1);
for j=1:nobs
    X{j} = [eye(nregions), repmat(groupRT(j),nregions,1)];
end

[BETA,sig,resid,vars,loglik2]  = mvregress(X,y2);
cc = [0.5 0.8 1; 1 0.6 0.8; 0.8 1 0.8];

B2 = [BETA(1:nregions)';repmat(BETA(end),1,nregions)];
xx = linspace(0.3,1.3)';
clf
for k = 1:3; 
plot(groupRT,groupEEG(:,k),'x','color',cc(k,:)); hold on
end

h = plot(xx, [ones(size(xx)),xx]*B2,'-','LineWidth',2);
legend(h,{sprintf('C1 beta: %s',num2str(BETA(1))),sprintf('C2 beta: %s',num2str(BETA(2))),sprintf('C3 beta: %s',num2str(BETA(3)))},'FontSize',8)
xlim([0 1.5]); ylim([-20 20])
title('EEG data at 3 sig clusters regressed against RT')
xlabel('RT'); ylabel('EEG')





