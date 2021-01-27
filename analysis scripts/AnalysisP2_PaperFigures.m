% Make the figures for the paper/thesis chapter
clear all; close all; clc; 
cd('W:\CIRCLE_EXP\LOG\EEG'); 
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1])
addpath('\\analyse2\Project0109\Lab\ckmatlab\eegck')
addpath('Z:\Lab\ckmatlab\ckinfo')
addpath('Z:\Lab\ckmatlab')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
load('GROUP_BEH_CUTOFF_19.mat'); 
n = 19; 

cc = [0.8 0.8 0.8];
Group_DIFF(:,1) = Group_RT(:,2) - Group_RT(:,1); 
Group_DIFF(:,2) = Group_RT(:,4) - Group_RT(:,3); 

Group_PerfD(:,1) = Group_Perf(:,2) - Group_Perf(:,1); 
Group_PerfD(:,2) = Group_Perf(:,4) - Group_Perf(:,3); 

subplot 121
stagger = 0.95:0.005:1.05; 
for k = 1:n; 
    plot(stagger(k),Group_DIFF(k,1),'o','MarkerFaceColor',[0.7 0.7 0.7]); hold on
    plot(stagger(k)+1,Group_DIFF(k,2),'o','MarkerFaceColor',[0.7 0.7 0.7]); hold on
end
hold on; boxplot(Group_DIFF); 
set(gca,'XTick',1:2,'XTickLabel',{'Auditory','Visual'})
ylabel('RT difference IC- C'); ylim([-0.3 0.3])
hline(0,':k'); 
title('Reaction Time Difference')


subplot 122
for k = 1:n; 
    plot(stagger(k),Group_PerfD(k,1),'o','MarkerFaceColor',[0.7 0.7 0.7]); hold on
    plot(stagger(k)+1,Group_PerfD(k,2),'o','MarkerFaceColor',[0.7 0.7 0.7]); hold on
end
hold on; boxplot(Group_PerfD); 
set(gca,'XTick',1:2,'XTickLabel',{'Auditory','Visual'})
ylabel('Performance difference IC- C'); 
ylim([-0.3 0.3])
hline(0,':k'); 
title('Performance Difference')

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_1_Behavioural_19.fig')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 : DECODER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding/GROUP_DECODE_CUTOFF_19.mat')
load('W:/CIRCLE_EXP/LOG/EEG/Files_AzBoot/Group_Az_19','prcntAUD','prcntVIS'); 

%--------------------------------------------------------------------------
% PLOT THE AZ
%--------------------------------------------------------------------------
TIMES = group_ND.TAX; 
clf
plot(TIMES,sq(mean(group_ND.audAz)),'LineWidth',3); hold on
plot(TIMES,sq(mean(group_ND.visAZ)),'r','LineWidth',3); 
legend('AUD','VIS')
vline(0,':k'); vline(0.3,':k')
ylabel('Az value'); xlabel('Time')
xlim([TIMES(1) TIMES(end)]);
title('Group Az')
hline([prcntAUD(2) prcntVIS(2)],{':b',':r'}); hold on
TAX = group_ND.TAX;

% add in where the Az stops being significant
audAz = sq(mean(group_ND.audAz));       % auditory Az
visAz = sq(mean(group_ND.visAZ));       % visual Az
ind = find(audAz>prcntAUD(2));          % ind AUD
indV = find(visAz>prcntVIS(2));         % ind VIS
vline(TAX(ind(113)),':b')                % auditory line
vline(TAX(indV(163)),':r')               % visual line

saveas(gcf,'Figure_3_Az_19.fig')


%--------------------------------------------------------------------------
% DECODER TOPOS
%--------------------------------------------------------------------------
load('W:/CIRCLE_EXP/LOG/EEG/Files_Decoding/GROUP_DECODE_CUTOFF_19.mat'); 

clf
% Plot Auditory Topo % Topoplot 
audTopo = mean(group_ND.audTopo,3);
visTopo = mean(group_ND.visTopo,3); 

subplot 121
load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat')                                          % load dummy file 
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
cfg.zlim = [-1 1];
% [~,J2] = max(sq(mean(group_ND.audAz)));                                      % plot the topography underlying the peak Az time point 
% EvpDummy.avg = sq(mean(mean(group_ND.audTopo(:,J2-2:J2+2,:),3),2));
ind = find(TIMES>=0.095 & TIMES<=0.2);
EvpDummy.avg = sq(mean(audTopo(:,ind),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('Auditory High vs Low \n %s to %s',num2str(TIMES(ind(1))),num2str(TIMES(ind(end))))); 

subplot 122
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.zlim = [-1 1];
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
% [~,J2] = max(sq(mean(group_ND.visAZ)));                                      % plot the topography underlying the peak Az time point 
% EvpDummy.avg = sq(mean(mean(group_ND.visTopo(:,J2-2:J2+2,:),3),2));
ind = find(TIMES>=0.145 & TIMES<=0.29);
EvpDummy.avg = sq(mean(visTopo(:,ind),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('Visual Small vs. Large \n %s to %s',num2str(TIMES(ind(1))),num2str(TIMES(ind(end))))); 

% suptitle('Group A projections, Max Az Points')

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_3_Decoding_Topos_19.fig')



%--------------------------------------------------------------------------
% DECODER FIPPED Y CONG VS INCONG WITH MC (AUD AND VIS)
%--------------------------------------------------------------------------
cd('W:\CIRCLE_EXP\LOG\EEG\Files_Decoding'); 
clear all; close all; clc;
load('Group_Shuffled_Ysignals_19'); 
load('Group_Shuffled_Ysignals_VIS_19'); 

addpath('W:\CIRCLE_EXP\LOG\EEG\functions')

% Plot it 
clf;
% suptitle('Group Y signals, corrected for multiple comparisons')
subplot 121; hold on;
plot(TAX,audCongruentM,'b','LineWidth',2); 
plot(TAX,audIncongruentM,'r','LineWidth',2); 
ind1 = find(audPC_Y.maskSig~=0); 
plot(TAX(ind1),audCongruentM(ind1),'k*'); 
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
hline(0,':k'); 
ylim([0.1 0.8]); xlabel('Time')
vline(0,':k'); vline(0.3,':k');
xlim([TAX(1) TAX(end)])
title('Auditory'); 


subplot 122; hold on;
plot(TAX,visCongruentM,'b','LineWidth',2); 
plot(TAX,visIncongruentM,'r','LineWidth',2); 
ind1 = find(visPC_Y.maskSig~=0); 
ind2 = find(visNC_Y.maskSig~=0); 
plot(TAX(ind1),visCongruentM(ind1),'k*'); 
plot(TAX(ind2),visCongruentM(ind2),'k*'); 
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
% hline(0,':k'); 
ylim([0.1 0.8])
vline(0,':k'); vline(0.3,':k'); 
title('Visual'); xlabel('Time')
xlim([TAX(1) TAX(end)])

suptitle('Group Y signals')
saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_3_DecodingSignals_MC.fig')



%--------------------------------------------------------------------------
% Regression with MCs
%--------------------------------------------------------------------------
clear all; close all; clc; 
load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Shuffled_RegWeights_19.mat')
load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Regression_19.mat','group_AW','group_AWD')

load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Shuffled_RegWeights_VIS_19.mat','TtrueV','Ttrue_VC')
load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Regression_VIS_19.mat','group_VW','group_VWD')

% load('Group_Regression_19.mat','TAX')
% I = find(TAX>=0 & TAX<=0.545);

% Take the regression weights from the smaller time window and plot them
% with the significant time points added in:
groupA = mean(group_AW(:,I)); 
groupC  = mean(group_AWD(:,I,1)); 
groupIC = mean(group_AWD(:,I,2)); 

groupV = mean(group_VW(:,I)); 
groupVC  = mean(group_VWD(:,I,1)); 
groupVIC = mean(group_VWD(:,I,2)); 

% TAX = TAX(I); 
suptitle('Group Regression Weights')

% One model, auditory neural weights 
clf; subplot 121; hold on
H1 = plot(TAX, groupA); 
H2 = plot(TAX, groupV,'r');                                                    % plot neural weights% plot neural weights
plot(TAX(oneModel.ind),groupA(oneModel.ind),'k*');              % plot significant points
title('One Model Regression Weights'); 
ylabel('Neural Weight'); legend([H1,H2], {'AUD','VIS'})
hline(0,':k'); vline(0,':k'); vline(0.3,':k'); 
xlabel('Time')
xlim([TAX(1) TAX(end)])

% One model, auditory tvalues 
subplot 122; hold on
H3 = plot(TAX,Ttrue,'color','b');                                % plot tvalues 
H4 = plot(TAX,TtrueV,'color','r');                                % plot tvalues 
h2 = plot(TAX(oneModel.ind),Ttrue(oneModel.ind),'k*');               % plot significant points 
title('One Model Regression TVALS'); 
xlabel('Time')
ylabel('t VALUE'); legend([H3,H4],{'AUD','VIS'})
hline(0,':k'); vline(0,':k'); vline(0.3,':k'); 
xlim([TAX(1) TAX(end)])

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_3_RegressionWeights_MC.fig')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MI FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc
cd('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient'); 
load('groupMI_gradient_Sig_VIS.mat')
load('W:\CIRCLE_EXP\LOG\EEG\128labels.mat'); 
labels = labels(1:128); 

subplot 221
imagesc(TAX(I),[],sq(groupAVG(:,:,1)),[0 0.1])
set(gca,'YTick',[32,64,96,128]); 
title('VIS congruent'); colorbar

subplot 222
imagesc(TAX(I),[],sq(groupAVG(:,:,2)),[0 0.1])
set(gca,'YTick',[32,64,96,128]); 
title('VIS incongruent'); colorbar

load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat'); 
load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/MI_gradient_stats_VIS.mat','vPos','vNeg'); 
TAX2 = TAX(I); 
subplot 223                                       % load dummy file 
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
cfg.zlim = [0 0.065];
ind = find(TAX2>=0.095 & TAX2<=0.2);
EvpDummy.avg = sq(mean(groupAVG(:,ind,1),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('VIS CONG \n %s to %s',num2str(TAX2(ind(1))),num2str(TAX2(ind(end))))); 

subplot 224                                      % load dummy file 
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
cfg.zlim = [0 0.05];
ind = find(TAX2>=0.095 & TAX2<=0.2);
EvpDummy.avg = sq(mean(groupAVG(:,ind,2),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('VIS INCONG \n %s to %s',num2str(TAX2(ind(1))),num2str(TAX2(ind(end))))); 

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_MI.fig')

%% compare congruent adn incongruent signals and get the regression weights 

clf
load('MI_gradient_stats_new.mat')
groupDiff = sq(mean(groupMI(:,:,:,1)-groupMI(:,:,:,2),3));  % congruent - incongruent

subplot(2,1,1)
plot(TAX(I), mean(groupAVG(:,:,1)),'LineWidth',2); hold on;
plot(TAX(I), mean(groupAVG(:,:,2)),'r','LineWidth',2);
hold on; 
ind = cat(1,aPos.clusterTimeInd{1},aPos.clusterTimeInd{2},aNeg.clusterTimeInd{1});
plot(TAX2(ind),mean(groupAVG(:,ind,1)),'k*')
xlim([0 0.545])
legend('Congruent','Incongruent')

EVP = EvpDummy; 
cfg = []; 
cfg.layout = 'biosemi128.lay'; 
cfg.comment = ' '; 
cfg.colorbar = 'yes'; 
cfg.zlim = [-0.02 0.02];
cfg.highlight = 'on';
cfg.highlightsymbol = '*'; 
cfg.highlightcolor = [0 0 0];

% np = 4;
% for k = 1:length(aPos.ind);
%     subplot(2,3,np);
%     EVP.avg = mean(groupDiff(:,aPos.clusterTimeInd{k}),2);
%     EVP.time = aPos.clusterTime{k}(1);
%     cfg.highlightchannel = aPos.clusterElecs{k};
%     ft_topoplotER(cfg,EVP)
%     title(sprintf('Auditory Positive Cluster \n %d, %s to %s',k,num2str(aPos.clusterTime{k}(1)),num2str(aPos.clusterTime{k}(end))));
% np = np+1;
% end

subplot 223
posInd = [aPos.clusterTimeInd{1} ;aPos.clusterTimeInd{2}];
posInd = unique(posInd); 
posElect = [aPos.clusterElecs{1} ;aPos.clusterElecs{2}];
EVP.avg = mean(groupDiff(:,posInd),2);
EVP.time = TAX2(posInd(1));
cfg.highlightchannel = posElect;
    ft_topoplotER(cfg,EVP)
title(sprintf('Auditory Positive Cluster \n %s to %s',num2str(TAX2(posInd(1))),num2str(TAX2(posInd(end)))));

% cfg.zlim = [-0.04 0.04];

subplot(2,2,4);
cfg.zlim = [-0.015 0.015];
EVP.avg = mean(groupDiff(:,aNeg.clusterTimeInd{1}),2);
% EVP.avg = mean(groupMI(:,aNeg.clusterTimeInd{1}),2);
EVP.time = aNeg.clusterTime{1}(1);
EVP.time = aNeg.clusterTime{1}(1);
cfg.highlightchannel = aNeg.clusterElecs{1};
ft_topoplotER(cfg,EVP)
title(sprintf('Negative Cluster \n %d, %s to %s',k,num2str(aNeg.clusterTime{1}(1)),num2str(aNeg.clusterTime{1}(end))));

% suptitle('Auditory Congruent MI - Auditory Incongruent MI')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_MI_comparisons.fig')


%% REGRESSION FIGURE
clearvars -except aPos aNeg
load('group_Regression_New')

clf
for s = 1:3;
    for k = 1:19;
        plot(s,groupB(k,s),'o'); hold on
    end
end
hold on; boxplot(groupB(:,1:3),'Colors','k');
ylabel('Beta Weight');
set(gca,'XTick',1:3,'XTickLabel',{'Cluster 1','Cluster 2','Cluster 3'})
title('Group Reg Weights: EEG:RT')

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_MI_Regression.fig')


%% VISUAL


clear all; clc
cd('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient'); 
load('groupMI_gradient_Sig_VIS.mat')
load('W:\CIRCLE_EXP\LOG\EEG\128labels.mat'); 
labels = labels(1:128); 

clf
set(gcf,'color','w')
subplot 211; hold on
plot(TAX(I),sq(mean(groupAVG(:,:,1),1))); 
plot(TAX(I),sq(mean(groupAVG(:,:,2),1)),'r'); 
xlim([TAX(I(1)) TAX(I(end))])
ylim([0 0.035])
legend('Congruent','Incongruent')
vline(0.07); vline(0.545); vline(0.06,'Color','r')

load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat'); 
load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/MI_gradient_stats_VIS.mat','vPos','vNeg'); 
TAX2 = TAX(I); 
subplot 223                                       % load dummy file 
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
cfg.zlim = [0 0.045];
ind = find(TAX2>=0.25 & TAX2<=0.35);
EvpDummy.avg = sq(mean(groupAVG(:,ind,1),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('VIS CONG \n %s to %s',num2str(TAX2(ind(1))),num2str(TAX2(ind(end))))); 

subplot 224                                      % load dummy file 
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.commentpos = 'title';
cfg.colorbar = 'yes';
% cfg.zlim = [0 0.05];
% ind = find(TAX2>=0.095 & TAX2<=0.2);
EvpDummy.avg = sq(mean(groupAVG(:,ind,2),2));
ft_topoplotER(cfg,EvpDummy)
title(sprintf('VIS INCONG \n %s to %s',num2str(TAX2(ind(1))),num2str(TAX2(ind(end))))); 



%%
clf
load('MI_gradient_stats_VIS.mat')
groupDiff = sq(mean(groupMI(:,:,:,1)-groupMI(:,:,:,2),3));  % congruent - incongruent

subplot(2,1,1)
plot(TAX(I), mean(groupAVG(:,:,1)),'LineWidth',2); hold on;
plot(TAX(I), mean(groupAVG(:,:,2)),'r','LineWidth',2);
hold on; 
ind = cat(1,vPos.clusterTimeInd{1},vNeg.clusterTimeInd{1},vNeg.clusterTimeInd{2});
plot(TAX2(ind),mean(groupAVG(:,ind,1)),'k*')
xlim([0 0.545])
legend('Congruent','Incongruent')


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
subplot 231
EVP.avg = mean(groupDiff(:,vPos.clusterTimeInd{1}),2);
EVP.time = vPos.clusterTime{1}(1);
cfg.highlightchannel = vPos.clusterElecs{1};
ft_topoplotER(cfg,EVP)
title(sprintf('Visual Positive Cluster %d, %s to %s',1,num2str(vPos.clusterTime{1}(1)),num2str(vPos.clusterTime{1}(end))));

for k = 1:length(vNeg.ind);
    subplot(2,3,k+1);
    EVP.avg = mean(groupDiff(:,vNeg.clusterTimeInd{k}),2);
    EVP.time = vNeg.clusterTime{k}(1);
    cfg.highlightchannel = vNeg.clusterElecs{k};
    ft_topoplotER(cfg,EVP)
    title(sprintf('Visual Negative Cluster %d, %s to %s',k,num2str(vNeg.clusterTime{k}(1)),num2str(vNeg.clusterTime{k}(end))));
end


load('group_Regression_New_VIS'); 
clf
for s = 1:3;
    for k = 1:19;
        plot(s,groupB(k,s),'o'); hold on
    end
end
hold on; boxplot(groupB(:,1:3),'Colors','k');
ylabel('Beta Weight');
set(gca,'XTick',1:3,'XTickLabel',{'Cluster 1','Cluster 2','Cluster 3'})
title('Group Reg Weights: EEG:RT VISUAL')




%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % MI Topos for significant time points where there was a C vc IC effect
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient'); 
% load 'groupMI_congruency'
% TAX2 = TAX(I); 
% 
% aTime{1} = find(TAX2>=0.08 & TAX2<=0.105);
% aTime{2} = find(TAX2>=0.14 & TAX2<=0.187);
% aTime{3} = find(TAX2>=0.209 & TAX2<=0.23);
% aTime{4} = find(TAX2>=0.312 & TAX2<=0.35);
% 
% groupAVG = mean(groupMI,4); 
% 
% load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat')
% clf
% subplot 221
% EVP = EvpDummy; 
% for k = 1:4; 
% EVP.avg = sq(mean(groupAVG(:,aTime{k},1),2));
% EVP.time = TAX(aTime{k}(1));
% 
% cfg = []; 
% cfg.zlim = [0 0.04];
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.colorbar = 'yes';
% title(sprintf('Auditory Cong %s to %s ms',num2str(TAX(aTime{k}(1))),num2str(TAX(aTime{k}(end))))) 
% ft_topoplotER(cfg,EVP)
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Mutual Information PLOTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all; clc; 
% load('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_congruency/groupMI_congruency.mat'); 
% 
% groupAVG = mean(groupMI,4); 
% maskAC = groupAVG(:,:,1).*IsigAC; 
% maskAIC = groupAVG(:,:,2).*IsigAIC; 
% maskVC = groupAVG(:,:,3).*IsigVC; 
% maskVIC = groupAVG(:,:,4).*IsigVIC; 
% 
% % MI/stimulus type figure
% subplot 221; imagesc(TAX(I),[],mean(groupMI(:,:,1,:),4),[0 0.06]); title('AUDITORY CONGRUENT'); colorbar
% subplot 222; imagesc(TAX(I),[],mean(groupMI(:,:,2,:),4),[0 0.06]); title('AUDITORY INCONGRUENT'); colorbar;
% subplot 223; imagesc(TAX(I),[],mean(groupMI(:,:,3,:),4),[0 0.06]); title('VISUAL CONGRUENT'); colorbar
% subplot 224; imagesc(TAX(I),[],mean(groupMI(:,:,4,:),4),[0 0.06]); title('VISUAL INCONGRUENT'); colorbar
% suptitle('Group Congruency MI')
% 
% % significance figure
% subplot 221; imagesc(TAX(I),[],IsigAC); title('AUDITORY CONGRUENT')
% subplot 222; imagesc(TAX(I),[],IsigAIC); title('AUDITORY INCONGRUENT')
% subplot 223; imagesc(TAX(I),[],IsigVC); title('VISUAL CONGRUENT')
% subplot 224; imagesc(TAX(I),[],IsigVIC); title('VISUAL INCONGRUENT')
% suptitle('Group Congruency MI significance')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %%
% % topoplots
% clf
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.zlim = [0 0.04];
% 
% %--------------------------------------------------------------------------
% % figure AC
% subplot 121; 
% ind = find(TAX2>=0.07 & TAX2<=0.135);
% EvpDummy.avg = mean(groupAVG(:,ind,1),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 122; 
% ind = find(TAX2>=0.15 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,1),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('AC'); colormap(redblue)
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % figure AIC
% clf
% subplot 121; 
% ind = find(TAX2>=0.07 & TAX2<=0.135);
% EvpDummy.avg = mean(groupAVG(:,ind,2),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 122; 
% ind = find(TAX2>=0.15 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,2),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('AIC'); colormap(redblue)
% %--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
% % figure VC
% clf
% cfg.zlim = [0 0.025];
% subplot 131; 
% ind = find(TAX2>=0.015 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 132; 
% ind = find(TAX2>=0.21 & TAX2<=0.250);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('VC'); colormap(redblue)
% 
% subplot 133; 
% ind = find(TAX2>=0.28 & TAX2<=0.32);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('VC'); colormap(redblue)
% %--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
% % figure VIC
% clf
% cfg.zlim = [0 0.025];
% subplot 131; 
% ind = find(TAX2>=0.015 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,4),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 132; 
% ind = find(TAX2>=0.21 & TAX2<=0.250);
% EvpDummy.avg = mean(groupAVG(:,ind,4),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('VC'); colormap(redblue)
% 
% subplot 133; 
% ind = find(TAX2>=0.28 & TAX2<=0.32);
% EvpDummy.avg = mean(groupAVG(:,ind,4),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('VIC'); colormap(redblue)
% %--------------------------------------------------------------------------
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% GROUP MEDIAN PLOTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clf
% clear all; clc;
% load('W:/CIRCLE_EXP/LOG/EEG/Files_MutualInformation/Files_MI_Median/groupMI_NaNs')
% order = [1 3 2 4]; 
% names = {'Auditory Stimulus','Visual Stimulus','Auditory Congruent','Visual Congruent'};
% for s = 1:4;
%     subplot(2,2,s); 
%     imagesc(TAX(I),[],g(:,:,order(s))); 
%     title(names{s})
% end
% saveas(gcf,'Figure_MedianMI.fig')
% %--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
% % significance plots 
% subplot 221; imagesc(TAX(I),[],IsigTone)
% title(names{1}); 
% 
% subplot 222; imagesc(TAX(I),[],IsigCircle)
% title(names{2}); 
% 
% subplot 223; imagesc(TAX(I),[],IsigAC)
% title(names{3})
% 
% subplot 224; imagesc(TAX(I),[],IsigVC); title('Visual Congruent')
% 
% saveas(gcf,'Figure_Median_Sig.fig')
% %--------------------------------------------------------------------------
% 
% load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat'); 
% clf
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.zlim = [0 0.04];
% 
% groupAVG = g; 
% TAX2 = TAX(I); 
% 
% %--------------------------------------------------------------------------
% % figure AC
% subplot 121; 
% ind = find(TAX2>=0.07 & TAX2<=0.135);
% EvpDummy.avg = mean(groupAVG(:,ind,1),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 122; 
% ind = find(TAX2>=0.15 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,1),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('AC'); colormap(redblue)
% saveas(gcf,'Figure_Topoplot_Aud.fig')
% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% % figure AIC
% clf
% cfg.zlim = [0 0.012];
% subplot 121; 
% ind = find(TAX2>=0.08 & TAX2<=0.11);
% EvpDummy.avg = mean(groupAVG(:,ind,2),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 122; 
% ind = find(TAX2>=0.34 & TAX2<=0.44);
% EvpDummy.avg = mean(groupAVG(:,ind,2),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('AIC'); colormap(redblue)
% saveas(gcf,'Figure_Topoplot_AC.fig')
% %--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
% % figure VC
% clf
% cfg.zlim = [0 0.025];
% subplot 131; 
% ind = find(TAX2>=0.015 & TAX2<=0.2);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% 
% subplot 132; 
% ind = find(TAX2>=0.21 & TAX2<=0.250);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('VC'); colormap(redblue)
% 
% subplot 133; 
% ind = find(TAX2>=0.28 & TAX2<=0.32);
% EvpDummy.avg = mean(groupAVG(:,ind,3),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy);
% colorbar
% suptitle('Visual Stim Median'); colormap(redblue)
% saveas(gcf,'Figure_Topoplot_Vis.fig')
% %--------------------------------------------------------------------------
% 
% 
% %--------------------------------------------------------------------------
% % figure VIC
% clf
% cfg.zlim = [0 0.01];
% % subplot 121; 
% ind = find(TAX2>=0.03 & TAX2<=0.36);
% EvpDummy.avg = mean(groupAVG(:,ind,4),2); 
% EvpDummy.time = mean(TAX2(ind)); 
% title(sprintf('[%s : %s]',num2str(TAX2(ind(1))),num2str(TAX2(ind(end)))))
% ft_topoplotER(cfg,EvpDummy); colorbar
% suptitle('Visual Congruent'); 
% saveas(gcf,'Figure_Topoplot_VC'); 
% 
% 
% 

%%
% 
% for k = 1:length(ind);
% elect{k} = EvpDummy.label(find(maskAC(:,ind(k))~=0));
% end
% 
% 
% %%
% % Find the significant time point/electrodes
% cfg = []; 
% cfg.zlim = [0 0.03];
% cfg.layout = 'biosemi128.lay';
% cfg.commentpos = 'title';
% cfg.comment = 'xlim';
% cfg.highlightchannel = elect;
% cfg.highlightsymbol    = '*';
% cfg.highlightcolor = [1 1 1];
% cfg.colorbar = 'yes';
% cfg.highlight = 'on';
% 
% % all significant time points AC
% clf
% [elec,time] = find(maskAC~=0);
% EvpDummy.avg = groupAVG(:,unique(time),1); 
% EvpDummy.time = TAX2(unique(time)); 
% cfg.xlim = EvpDummy.time(3):0.02:EvpDummy.time(end);
% ft_topoplotER(cfg,EvpDummy);
% suptitle('AC')
% 
% clf
% % all significant time points AIC
% [elec,time] = find(maskAIC~=0);
% EvpDummy.avg = groupAVG(:,unique(time),2); 
% EvpDummy.time = TAX2(unique(time)); 
% cfg = []; 
% cfg.layout = 'biosemi128.lay';
% cfg.commentpos = 'title';
% cfg.comment = 'xlim';
% cfg.highlightchannel = chans;
% cfg.highlightsymbol    = '+';
% cfg.highlightcolor = [1 1 1];
% cfg.colorbar = 'yes';
% cfg.xlim = EvpDummy.time(3):0.02:EvpDummy.time(end);
% ft_topoplotER(cfg,EvpDummy);
% suptitle('AIC')
% 
% clf
% % all significant time points VC
% [elec,time] = find(maskVC~=0);
% EvpDummy.avg = groupAVG(:,unique(time),3); 
% EvpDummy.time = TAX2(unique(time)); 
% cfg = []; 
% cfg.layout = 'biosemi128.lay';
% cfg.commentpos = 'title';
% cfg.comment = 'xlim';
% cfg.highlightchannel = chans;
% cfg.highlightsymbol    = '+';
% cfg.highlightcolor = [1 1 1];
% cfg.colorbar = 'yes';
% cfg.xlim = EvpDummy.time(3):0.04:EvpDummy.time(end);
% ft_topoplotER(cfg,EvpDummy);
% suptitle('VC')
% 
% 
% clf
% % all significant time points VIC
% [elec,time] = find(maskVIC~=0);
% EvpDummy.avg = maskVIC(:,unique(time)); 
% EvpDummy.time = TAX2(unique(time)); 
% cfg = []; 
% cfg.layout = 'biosemi128.lay';
% cfg.commentpos = 'title';
% cfg.comment = 'xlim';
% cfg.highlightchannel = chans;
% cfg.highlightsymbol    = '+';
% cfg.highlightcolor = [1 1 1];
% cfg.colorbar = 'yes';
% cfg.xlim = EvpDummy.time(3):0.04:EvpDummy.time(end);
% ft_topoplotER(cfg,EvpDummy);
% suptitle('VIC')
% 
% 
% 
% %%
% 
% 




























%% SUPPLEMENTARY FIGURES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure X. Decoder Separation Projections : Result, not the flip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; 

load('Group_Ysignals_19'); 

% High tone and Low tone Y signal Plots with individual subjects
clf; subplot 221; hold on
for k = 1:n; 
    plot(TIMES,mean(HighToneG{k}(:,I)),'Color',[1 0.7 0.7]); hold on
    plot(TIMES,mean(LowToneG{k}(:,I)),'Color',[0 0.8 1]); hold on
end
h1 = plot(TIMES,mean(meanHT(:,I)),'r','LineWidth',3);
h2 = plot(TIMES,mean(meanLT(:,I)),'b','LineWidth',3);
ylim([-3 3])
xlim([TIMES(1) TIMES(end)]); vline(0.3,':k');
legend([h1 h2],'HIGH TONE','LOW TONE'); 
ylabel('AUD Y SIGNALS'); title('NON FLIPPED');

% Congruent vs Incongruent Y signal Plots with individual subjects
subplot 222; hold on
for k = 1:n; 
    plot(TIMES,mean(CongruentG{k}(:,I)),'Color',[0 0.8 0.8]); hold on
    plot(TIMES,mean(IncongruentG{k}(:,I)),'Color',[1 0.5 0.8]); hold on
end
h1 = plot(TIMES,mean(meanC(:,I)),'g','LineWidth',3);
h2 = plot(TIMES,mean(meanIC(:,I)),'m','LineWidth',3);
ylim([0-0.35 0.35])
xlim([TIMES(1) TIMES(end)]); vline(0.3,':k');
legend([h1 h2],'CONG','INCONG')
ylabel('AUD Y SIGNALS');title('NON FLIPPED');

% VISUAL Y SIGNALS 
load('Group_Ysignals_VIS_19')
I = find(TIMES>=0 & TIMES<=0.545); 
TIMES_save = TIMES; 
TIMES = TIMES(I);
subplot 223; hold on
for k = 1:n; 
    plot(TIMES,mean(SmallCircleG{k}(:,I)),'Color',[1 0.7 0.7]); hold on
    plot(TIMES,mean(LargeCircleG{k}(:,I)),'Color',[0 0.8 1]); hold on
end
h1 = plot(TIMES,mean(meanSC(:,I)),'r','LineWidth',3);
h2 = plot(TIMES,mean(meanLC(:,I)),'b','LineWidth',3);
ylim([-3 3])
xlim([TIMES(1) TIMES(end)]); vline(0.3,':k');
legend([h1 h2],'Small Circle','Large Circle'); 
ylabel('VIS Y SIGNALS'); title('NON FLIPPED');

subplot 224; hold on
for k = 1:n; 
    plot(TIMES,mean(CongruentG{k}(:,I)),'Color',[0 0.8 0.8]); hold on
    plot(TIMES,mean(IncongruentG{k}(:,I)),'Color',[1 0.5 0.8]); hold on
end
h1 = plot(TIMES,mean(meanC(:,I)),'g','LineWidth',3);
h2 = plot(TIMES,mean(meanIC(:,I)),'m','LineWidth',3);
ylim([0-0.35 0.35])
xlim([TIMES(1) TIMES(end)]); vline(0.3,':k');
legend([h1 h2],'CONG','INCONG')
ylabel('VIS Y SIGNALS'); title('NON FLIPPED');

saveas(gcf,'Figure_4_Ysignals_19.fig')



% 
saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_Supp_Decoding_19.fig')


%%
clear all; close all; clc;
load('GROUP_DECODE_CUTOFF_19.mat')
load('W:/CIRCLE_EXP/LOG/EEG/Files_AzBoot/Group_Az_19','prcntAUD','prcntVIS'); 

TIMES = group_ND.TAX; 
clf
subplot 121
for s = 1:19; 
    plot(TIMES,group_ND.audAz(s,:),'Color',[0 0.8 1]); hold on
end
plot(TIMES,sq(mean(group_ND.audAz)),'LineWidth',3); hold on
legend('AUD')
vline(0,':k'); vline(0.3,':k')
ylabel('Az value'); xlabel('Time')
xlim([TIMES(1) TIMES(end)]);
title('Group Az')
hline([prcntAUD(2)],':b'); hold on
audAz = sq(mean(group_ND.audAz));       % auditory Az
ind = find(audAz>prcntAUD(2));          % ind AUD
vline(TAX(ind(113)),':b')                % auditory line

subplot 122
for s = 1:19
        plot(TIMES,group_ND.visAZ(s,:),'Color',[1 0.7 0.7]); hold on;
end
plot(TIMES,sq(mean(group_ND.visAZ)),'r','LineWidth',3); 
legend('VIS')
vline(0,':k'); vline(0.3,':k')
ylabel('Az value'); xlabel('Time')
xlim([TIMES(1) TIMES(end)]);
title('Group Az')
hline([prcntVIS(2)],':r'); hold on

TAX = group_ND.TAX;

% add in where the Az stops being significant
visAz = sq(mean(group_ND.visAZ));       % visual Az
indV = find(visAz>prcntVIS(2));         % ind VIS
vline(TAX(indV(163)),':r')               % visual line

saveas(gcf,'Figure_Supplementary_Fig1_Az.fig')



%%

cd('W:\CIRCLE_EXP\LOG\EEG'); 
clear all; close all; clc;
load('Group_Shuffled_Ysignals_19'); 
load('Group_Shuffled_Ysignals_VIS_19'); 

addpath('W:\CIRCLE_EXP\LOG\EEG\functions')

% Plot it 
clf;
% suptitle('Group Y signals, corrected for multiple comparisons')
subplot 121; hold on;

for s = 1:19;
    plot(TAX,audCongruent(s,:),'color',[0 0.8 1]); hold on

end
for s = 1:19;
    plot(TAX,audIncongruent(s,:),'color',[1 0.7 0.7]); hold on
end
plot(TAX,audCongruentM,'b','LineWidth',3); 
plot(TAX,audIncongruentM,'r','LineWidth',3); 

ind1 = find(audPC_Y.maskSig~=0); 
plot(TAX(ind1),audCongruentM(ind1),'k*'); 
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
hline(0,':k'); 
ylim([0.1 0.8]); xlabel('Time')
vline(0,':k'); vline(0.3,':k');
xlim([TAX(1) TAX(end)])
title('Auditory'); 


subplot 122; hold on;
plot(TAX,visCongruentM,'b','LineWidth',2); 
plot(TAX,visIncongruentM,'r','LineWidth',2); 
ind1 = find(visPC_Y.maskSig~=0); 
ind2 = find(visNC_Y.maskSig~=0); 
plot(TAX(ind1),visCongruentM(ind1),'k*'); 
plot(TAX(ind2),visCongruentM(ind2),'k*'); 
legend('Congruent','Incongruent','p<0.05')
ylabel('Y FLIPPED')
% hline(0,':k'); 
ylim([0.1 0.8])
vline(0,':k'); vline(0.3,':k'); 
title('Visual'); xlabel('Time')
xlim([TAX(1) TAX(end)])

suptitle('Group Y signals')
saveas(gcf,'Figure_4_DecodingSignals_MC.fig')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regression with MCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; 

load('W:/CIRCLE_EXP/LOG/EEG/Group_Regression_19.mat','group_AW','group_AWD')
load('W:/CIRCLE_EXP/LOG/EEG/Group_Regression_VIS_19.mat','group_VW','group_VWD')
load('W:/CIRCLE_EXP/LOG/EEG/Group_Shuffled_Ysignals_VIS_19.mat','TAX','I')

clf
subplot 121; 
for s = 1:19; 
    plot(TAX,group_AW(s,I),'Color',[0.7 0.7 0.7]); hold on
end
groupA = mean(group_AW(:,I)); 
hold on; plot(TAX,groupA,'LineWidth',3); 
load('Group_Shuffled_RegWeights_19.mat'); 
plot(TAX(oneModel.ind),groupA(oneModel.ind),'k*');     

ylim([-1.5 1.7]); xlim([TAX(1) TAX(end)]); ylabel('Auditory Neural Weights'); 
subplot 122; 
for s = 1:19; 
    plot(TAX,group_VW(s,I),'Color',[0.7 0.7 0.7]); hold on
end
hold on; plot(TAX,mean(group_VW(:,I)),'r','LineWidth',3); 
ylim([-1.5 1.7]); xlim([TAX(1) TAX(end)]); ylabel('Visual Neural Weights'); 

% load('Group_Regression_19.mat','TAX')
% I = find(TAX>=0 & TAX<=0.545);

% Take the regression weights from the smaller time window and plot them
% with the significant time points added in:
groupA = mean(group_AW(:,I)); 


groupV = mean(group_VW(:,I)); 


% TAX = TAX(I); 
suptitle('Group Regression Weights')

% One model, auditory neural weights 
clf; subplot 121; hold on
H1 = plot(TAX, groupA) 
H2 = plot(TAX, groupV,'r')                                                    % plot neural weights% plot neural weights
plot(TAX(oneModel.ind),groupA(oneModel.ind),'k*');              % plot significant points
title('One Model Regression Weights'); 
ylabel('Neural Weight'); legend([H1,H2], {'AUD','VIS'})
hline(0,':k'); vline(0,':k'); vline(0.3,':k'); 
xlabel('Time')
xlim([TAX(1) TAX(end)])

% One model, auditory tvalues 
subplot 122; hold on
H3 = plot(TAX,Ttrue,'color','b')                                % plot tvalues 
H4 = plot(TAX,TtrueV,'color','r')                                % plot tvalues 
h2 = plot(TAX(oneModel.ind),Ttrue(oneModel.ind),'k*');               % plot significant points 
title('One Model Regression TVALS'); 
xlabel('Time')
ylabel('t VALUE'); legend([H3,H4],{'AUD','VIS'})
hline(0,':k'); vline(0,':k'); vline(0.3,':k'); 
xlim([TAX(1) TAX(end)])




%% MUTUAL INFORMATION FIGURE

% part 1: MI values 
load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/groupMI')
load('W:\CIRCLE_EXP\LOG\EEG/EVP128.mat')

group_miAC_avg = mean(groupmiAC,3); 
group_miAIC_avg = mean(groupmiAIC,3); 
group_diff = group_miAC_avg - group_miAIC_avg; 

subplot 221; 
imagesc(TAX(I),[],sq(mean(groupmiAC,3))); 
ylabel('Electrode'); xlabel('Time (ms)'); 
title('CONGRUENT MI'); 

subplot 222; 
imagesc(TAX(I),[],sq(mean(groupmiAIC,3))); 
ylabel('Electrode'); xlabel('Time (ms)'); 
title('INCONGRUENT MI'); 

subplot 223; 
ind = find(TAX2>=0.1 & TAX2<=0.2); 
EvpDummy.avg = sq(mean(group_miAC_avg(:,ind),2)); 
cfg = []; 
cfg.layout = 'biosemi128.lay'; 
cfg.comment = ' ';
cfg.zlim = [0 0.06];
ft_topoplotER(cfg,EvpDummy)
colorbar
title('AUD CONG MI 0.1:0.2'); 

subplot 224; 
ind = find(TAX2>=0.1 & TAX2<=0.2); 
EvpDummy.avg = sq(mean(group_miAIC_avg(:,ind),2)); 
cfg = []; 
cfg.layout = 'biosemi128.lay'; 
cfg.comment = ' ';
cfg.zlim = [0 0.06];
ft_topoplotER(cfg,EvpDummy)
colorbar
title('AUD INCONG MI 0.1:0.2'); 

load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/tvalues'); 
load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/MI_gradient_stats_new.mat')
ind_all = [aPos.clusterTimeInd{1}; aPos.clusterTimeInd{2}; aNeg.clusterTimeInd{1}]; 
ind_all = unique(ind_all); 
TAX2 = TAX(I); 
plot(TAX2,mean(TtrueA),'LineWidth',2); hold on
plot(TAX2(ind_all),mean(TtrueA(:,ind_all)),'*r')
xlim([0 0.545]); ylim([-1.4 1.4]); 
hline(0)

saveas(gcf,'W:/CIRCLE_EXP/LOG/EEG/Figures_Thesis/Figure_4_MI.fig');

%%
% Figure 5.

load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/MI_gradient_stats_new.mat')
load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat');
groupDiff = sq(mean(groupMI(:,:,:,1)-groupMI(:,:,:,2),3));

EVP = EvpDummy;
cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.comment = ' ';
cfg.colorbar = 'yes';
cfg.zlim = [-0.02 0.02];
cfg.highlight = 'on';
cfg.highlightsymbol = '*';
cfg.highlightcolor = [0 0 0];
clf
np = 1; 
for k = 1:length(aPos.ind);
    subplot(2,3,np);
    EVP.avg = mean(groupDiff(:,aPos.clusterTimeInd{k}),2);
    EVP.time = aPos.clusterTime{k}(1);
    cfg.highlightchannel = aPos.clusterElecs{k};
    ft_topoplotER(cfg,EVP)
    title(sprintf('Auditory Positive Cluster %d, %s to %s',k,num2str(aPos.clusterTime{k}(1)),num2str(aPos.clusterTime{k}(end))));
    np = np+1; 
end

cfg.zlim = [-0.04 0.04];

subplot(2,3,np);
cfg.zlim = [-0.02 0.02];
EVP.avg = mean(groupDiff(:,aNeg.clusterTimeInd{1}),2);
EVP.time = aNeg.clusterTime{1}(1);
EVP.time = aNeg.clusterTime{1}(1);
cfg.highlightchannel = aNeg.clusterElecs{1};
ft_topoplotER(cfg,EVP)
title(sprintf('Negative Cluster %d, %s to %s',k,num2str(aNeg.clusterTime{1}(1)),num2str(aNeg.clusterTime{1}(end))));


load('W:\CIRCLE_EXP\LOG\EEG\Files_MutualInformation\Files_MI_Gradient/group_Regression_New'); 

subplot(2,1,2)
for k = 1:3; 
    for s = 1:19; 
        plot(k,groupB(s,k),'o'); hold on
    end
end
boxplot(groupB(:,1:3),'Colors','k')
xlim([0.5 3.5]);
set(gca,'XTick',1:3,'XTickLabel',{'Cluster 1','Cluster 2','Cluster 3'}); 
ylabel('Beta Weight'); 


saveas(gcf,'W:\CIRCLE_EXP\LOG\EEG\Figures_Thesis/Figure_5_SigMI_Regression_NEW.fig')
%%
% 
% %%
% indAud = [81 82]; 
% indVis = [61:63 106:108 119:120 166:170];
% 
% v{1} = 61:63; v{2} = 106:108; v{3} = 119:120; v{4} = 166:170; 
% 
% 
% %% 
% clear all; close all; clc 
% 
% AUD = []; VIS = []; aTrials = []; vTrials = []; 
% for s = 1:20; 
%    
%    load(sprintf('W:/CIRCLE_EXP/LOG/EEG/Analysis_MutualInformation/Files_DataPrepro/S0%d_EEG_matrix.mat',s));
%     if s==1; 
%         I = find(TAX>=0 & TAX<=0.545); 
%     end
%     
%    AUD = cat(1,AUD,data(:,:,I));
%    aTrials = cat(1,aTrials,trials);
%    
%    load(sprintf('W:/CIRCLE_EXP/LOG/EEG/Analysis_MutualInformation/S0%d_EEG_matrix_VIS.mat',s));
%    VIS = cat(1,VIS,data(:,:,I)); 
%    vTrials = cat(1,vTrials,trials); 
%    fprintf('S0%d \n',s); 
%    
% end
% %%
% cd('W:\CIRCLE_EXP\LOG\EEG'); 
% 
% AC = sq(mean(AUD(aTrials(:,3)==1,:,:)));
% AIC = sq(mean(AUD(aTrials(:,3)==2,:,:)));
% 
% VC = sq(mean(VIS(vTrials(:,3)==1,:,:)));
% VIC = sq(mean(VIS(vTrials(:,3)==2,:,:)));
% 
% 
% indAud = [21 22]; 
% v{1} = 1:3; v{2} = 46:48; v{3} = 59:60; v{4} = 106:110; 
% 
% load('EVP128.mat')   
% %%
% clf
% % Plot Auditory Topo % Topoplot Difference 
% 
% AudDiff = [sq(mean(AIC(:,indAud),2)) - sq(mean(AC(:,indAud),2))]; % IC - C
% 
% 
% subplot 151 % load dummy file 
% cfg = [];
% cfg.layout = 'biosemi128.lay';
% cfg.comment = ' ';
% cfg.commentpos = 'title';
% cfg.colorbar = 'yes';
% cfg.zlim = [-0.25 0.25];                                 % plot the topography underlying the peak Az time point 
% EvpDummy.avg = AudDiff;
% ft_topoplotER(cfg,EvpDummy); 
% title('AUD INCONGRUENT - CONGRUENT, 100:105 ms')
% 
% set(gcf,'color','w')
% saveas(gcf,'Figure_8_topoYaud_Diffs.fig')
% 
% %% VISUAL 
% TAX2 = TAX(I);
% for k = 1:4; 
% VisDiff(:,k) = [sq(mean(VIC(:,v{k}),2)) - sq(mean(VC(:,v{k}),2))]; % IC - C
% end
% 
% cfg = [];
% cfg.layout = 'biosemi128.lay';
% cfg.comment = ' ';
% cfg.commentpos = 'title';
% cfg.colorbar = 'yes';
% cfg.zlim = [-0.25 0.25];    
% 
% 
% for k = 1:4
% subplot(1,5,k+1)                               
% EvpDummy.avg = VisDiff(:,k);                  % plot the topography underlying the peak Az time point
% ft_topoplotER(cfg,EvpDummy)
% 
% title(sprintf('VIS IC-C, %s:%s ms',num2str(TAX2(v{k}(1))),num2str(TAX2(v{k}(end)))))
% end
% 
% saveas(gcf,'Figure_8_topoYvis_Diffs.fig')
% saveas(gcf,'Figure_9_topos_YDiffs.fig')





% 
% 
% %% 
% % TOPOS UNDELRYING THE SIGNIFICANT DECODING TIME POINTS 
% 
% clear all; close all; clc; 
% cd('W:/CIRCLE_EXP/LOG/EEG');
% addpath('Z:\Lab\ckmatlab\ckinfo')
% load('W:/CIRCLE_EXP/LOG/EEG/Files_Decoding/GROUP_DECODE_CUTOFF_19.mat')
% 
% I = find(TIMES>=0 & TIMES<=0.545); 
% audTopo = group_ND.audTopo(:,I,:);
% visTopo = group_ND.visTopo(:,I,:);
% 
% 
% TAX2 = TIMES(I); 
% acTimes = find(TAX2>=0.095 & TAX2<=0.110);
% vcTimes{1} = find(TAX2>=0.200 & TAX2<=0.235);
% vcTimes{2} = find(TAX2>=0.460 & TAX2<=0.465);
% vcTimes{3} = find(TAX2>=0.530 & TAX2<=0.545);
% 
% 
% 
% load('EVP128.mat'); 
% 
% % Auditory ones 
% EVP = EvpDummy; 
% EVP.avg = sq(mean(mean(audTopo(:,acTimes,:),3),2));
% EVP.time = mean(TAX2(acTimes));
% 
% 
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.colorbar = 'yes';
% ft_topoplotER(cfg,EVP)
% title(sprintf('AUD TOPO, CONG vs. INCONG Time point \n %s to %s',num2str(TAX2(acTimes(1))),num2str(TAX2(acTimes(end)))))
% saveas(gcf,'Figure_AudCong_Topo_SigTimepoints.fig')
% 
% 
% % Visual
% EVP = EvpDummy; 
% for k = 1:3; 
%     clf
% EVP.avg = sq(mean(mean(visTopo(:,vcTimes{k},:),3),2));
% EVP.time = mean(TAX2(acTimes));
% 
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.colorbar = 'yes';
% ft_topoplotER(cfg,EVP)
% title(sprintf('VIS TOPO, CONG vs. INCONG Time point \n %s to %s',num2str(TAX2(vcTimes{k}(1))),num2str(TAX2(vcTimes{k}(end)))))
% saveas(gcf,sprintf('Figure_VisCong%d_Topo_SigTimepoints.fig',k))
% end
% 
% 
% 
% %% REGRESSION SIGNIFICANT TIME POINTS
% clear acTimes
% load('W:/CIRCLE_EXP/LOG/EEG/Files_Decoding/GROUP_DECODE_CUTOFF_19.mat')
% 
% I = find(TIMES>=0 & TIMES<=0.545); 
% audTopo = group_ND.audTopo(:,I,:);
% visTopo = group_ND.visTopo(:,I,:);
% 
% acTimes{1} = find(TAX2>=0.065 & TAX2<=0.135);
% acTimes{2} = find(TAX2>=0.180 & TAX2<=0.205);
% 
% % first one
% EVP = EvpDummy; 
% EVP.avg = sq(mean(audTopo(:,acTimes{1},:),3));
% EVP.time = TAX2(acTimes{1});
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = 'xlim'; 
% cfg.commentpos = 'title';
% cfg.colorbar = 'yes';
% cfg.xlim = EVP.time; 
% ft_topoplotER(cfg,EVP)
% suptitle('AUD REGRESSION early time points Topos')
% saveas(gcf,'Figure_AudReg1_Topo_SigTimepoints.fig')
% 
% 
% EVP = EvpDummy; 
% clf
% EVP.avg = sq(mean(audTopo(:,acTimes{2},:),3));
% EVP.time = TAX2(acTimes{2});
% cfg = []; 
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = 'xlim'; 
% cfg.commentpos = 'title';
% cfg.colorbar = 'yes';
% cfg.xlim = EVP.time; 
% ft_topoplotER(cfg,EVP)
% suptitle('AUD REGRESSION later time points Topos')
% saveas(gcf,'Figure_AudReg2_Topo_SigTimepoints.fig')
% 
% 






% 
% %--------------------------------------------------------------------------
% % REG WEIGHTS AND TVALS MODALITY MODEL WITH MC : AUD AND VIS 
% %--------------------------------------------------------------------------
% addpath('W:\CIRCLE_EXP\LOG\EEG\functions')
% load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Regression_19.mat')
% load('W:/CIRCLE_EXP/LOG/EEG/Files_Regression/Group_Regression_VIS_19.mat')
% 
% I = find(Decode_times>=0 & Decode_times<=0.5);          % narrow down the time window (-0.3 to 0.5 stimulus onset)
% TAX = Decode_times(I);                                  % new timings
% 
% aw = group_AW;                                          % one model, auditory weight 
% awD = group_AWD;                                        % dual model, modelled separately by congruency;
% vw = group_VW;                                          % one model, auditory weight 
% vwD = group_VWD;                                        % dual model, modelled separately by congruency;
% 
% 
% % Plot : auditory regression weights one model
% clf
% subplot 231; hold on
% plot(TAX,mean(aw(:,I)),'LineWidth',2,'Color','b');
% plot(TAX,mean(vw(:,I)),'LineWidth',2,'Color','r'); 
% ylim([-0.3 0.3])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% title('Neural Weights: Single Model'); ylabel('REGRESSION WEIGHT'); 
% legend('AUD','VIS')
% xlabel('Time')
% 
% 
% % Plot: auditory regression weights dual model 
% subplot 232; hold on 
% plot(TAX,mean(awD(:,I,1)),'LineWidth',2,'Color','b');   % congruent
% plot(TAX,mean(awD(:,I,2)),'LineWidth',2,'Color','m')    % incongruent
% ylim([-0.3 0.3])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% title('A Weights: Dual Model'); ylabel('REGRESSION WEIGHT')
% legend('CONG','INCONG')
% xlabel('Time')
% 
% 
% subplot 233; hold on
% plot(TAX,mean(vwD(:,I,1)),'LineWidth',2,'Color','k');   % congruent
% plot(TAX,mean(vwD(:,I,2)),'LineWidth',2,'Color','r')    % incongruent
% ylim([-0.3 0.3])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% title('V Weights: Dual Model'); ylabel('REGRESSION WEIGHT')
% legend('CONG','INCONG')
% xlabel('Time')
% 
% 
% % Calculate T values for each
% tA = zeros(1,size(aw,2)); tAC = tA; tAIC = tA;
% tV = zeros(1,size(aw,2)); tVC = tV; tVIC = tV;
% for k = 1:length(aw);
%     tA(k) = mean(aw(:,k))./sem(aw(:,k));            % 1 model : auditory
%     tAC(k) = mean(awD(:,k,1))./sem(awD(:,k,1));     % Dual model : auditory congruent
%     tAIC(k) = mean(awD(:,k,1))./sem(awD(:,k,2));    % Dual model : auditory incongruent
%     
%     tV(k) = mean(vw(:,k))./sem(vw(:,k));            % 1 model : auditory
%     tVC(k) = mean(vwD(:,k,1))./sem(vwD(:,k,1));     % Dual model : auditory congruent
%     tVIC(k) = mean(vwD(:,k,1))./sem(vwD(:,k,2));    % Dual model : auditory incongruent
%     
% end
% 
% % One model : auditory t vals
% subplot 234; hold on
% shadedErrorBar(TAX,tA(I),sem(aw(:,I)),{'b','markerfacecolor','b'})
% shadedErrorBar(TAX,tV(I),sem(vw(:,I)),{'r','markerfacecolor','r'})
% ylim([-4 4])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% ylabel('TVALUE')
% title('TVALS, all Aud, bars = sem')
% xlabel('Time')
% 
% 
% % Dual model : auditory tvals
% subplot 235; hold on;
% shadedErrorBar(TAX,tAC(I),sem(awD(:,I,1)),{'b','markerfacecolor','b'})      % congruent
% shadedErrorBar(TAX,tAIC(I),sem(awD(:,I,2)),{'m','markerfacecolor','m'})     % incongruent 
% ylim([-4 4])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% ylabel('TVALUE')
% title('TVALS, all Aud, bars = sem')
% xlabel('Time')
% 
% 
% % Dual model : auditory tvals
% subplot 236; hold on;
% shadedErrorBar(TAX,tVC(I),sem(awD(:,I,1)),{'b','markerfacecolor','b'})      % congruent
% shadedErrorBar(TAX,tVIC(I),sem(awD(:,I,2)),{'m','markerfacecolor','m'})     % incongruent 
% ylim([-4 4])
% vline(0,':k'); vline(0.3,':k'); hline(0,':k')
% ylabel('TVALUE')
% title('TVALS, all VIS, bars = sem'); legend('CONG','INCONG')
% suptitle('Group regression weights'); xlabel('Time')
% 
% saveas(gcf,'Figure_3_RegressionWeights_19.fig')

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% FIGURE 7 : CORRELATIONS WITH MC : AUD AND VIS 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear all;
% close all; clc 
% load('W:/CIRCLE_EXP/LOG/EEG/Group_Correlations_Corrected_19.mat')
% 
% % clear TAX; TAX = Decode_times(I); 
% 
% clf; hold on;
% H1 = plot(TAX,rAUD,'b','LineWidth',2); 
% ind1 = find(PC.maskSig~=0); 
% ind2 = find(NC.maskSig~=0); 
% plot(TAX(ind1),rAUD(ind1),'k*'); 
% plot(TAX(ind2),rAUD(ind2),'k*'); 
% ylabel('R VALUE'); 
% legend(H1,{'AUD'});
% hline(0); vline(0.3); 
% suptitle('Correlation (AIC-C) RT with (AIC-C) Flipped Y signals')
% title('Non Bootstrapped Trials (non equal)')
% 
% 
% hold on
% load('W:/CIRCLE_EXP/LOG/EEG/Group_Correlations_Corrected_VIS_19.mat')
% H2 = plot(TAX,rVIS,'r','LineWidth',2); hold on
% ind1 = find(PC.maskSig~=0); 
% ind2 = find(NC.maskSig~=0); 
% plot(TAX(ind1),rVIS(ind1),'k*'); 
% plot(TAX(ind2),rVIS(ind2),'k*'); 
% ylabel('R VALUE'); 
% legend([H1,H2],{'AUD','VIS'});
% hline(0); vline(0.3); 
% suptitle('Correlation (AIC-C) RT with (AIC-C) Flipped Y signals')
% title('Non Bootstrapped Trials (non equal)')
% 
% xlim([TAX(1) TAX(end)])
% 
% saveas(gcf,'Figure_7_Correlations_19.fig')


% % add topo that goes with the significant time points 
% load('W:/CIRCLE_EXP/LOG/EEG/Files_Decoding/GROUP_DECODE_CUTOFF_19.mat','group_ND')
% 
% audTopo = group_ND.audTopo(:,I,:);
% visTopo = group_ND.visTopo(:,I,:);
% acTimes{1} = find(TAX>=0.065 & TAX<=0.135);
% acTimes{2} = find(TAX>=0.180 & TAX<=0.205);
% 
% % first one
% load('W:/CIRCLE_EXP/LOG/EEG/EVP128.mat')
% clf
% subplot 221
% EVP = EvpDummy; 
% EVP.avg = sq(mean(mean(audTopo(:,acTimes{1},:),3),2));
% EVP.time = TAX(acTimes{1}(1));
% cfg = []; 
% cfg.zlim = [-2 2];
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.colorbar = 'yes';
% title(sprintf('Auditory %s to %s ms',num2str(TAX(acTimes{1}(1))),num2str(TAX(acTimes{1}(end))))) 
% ft_topoplotER(cfg,EVP)
% 
% % later one 
% subplot 222 
% EVP = EvpDummy; 
% EVP.avg = sq(mean(mean(audTopo(:,acTimes{2},:),3),2));
% EVP.time = TAX(acTimes{2}(1));
% cfg = []; 
% cfg.zlim = [-2 2];
% cfg.layout = 'biosemi128.lay'; 
% cfg.comment = ' '; 
% cfg.colorbar = 'yes';
% title(sprintf('Auditory %s to %s ms',num2str(TAX(acTimes{2}(1))),num2str(TAX(acTimes{2}(end))))) 
% ft_topoplotER(cfg,EVP)
% 
% suptitle('Group Decoding A topos, underlying the significant time points of Aud reg weights')
% saveas(gcf,'Figure_DecodingTopos_SigRegression_TimePoints_AUD.fig')

