addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

% define the cohort
sbj_names_all = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
    ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
    ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
    ;'S17_116_AA';'S17_118_TW';'S19_145_PC';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};

%make a specific selection of cohort
indxcohort = 1:36;%China
% indxcohort = [37:44];%Stanfordc
% indxcohort = [1:44];%2 Centers

sbj_names = sbj_names_all(indxcohort);%China


% define the the abbreviations of kinds of brian structures
anat_all = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','WM','OUT','EC',...
    'FOP','POP','TOP','EMPTY','PARL','LESION','INSULA','BASAL'};

%make a specific selection of anatomical structures
anat = {'HIPPO A','HIPPO M','HIPPO P','PHG'};anat_name = 'MTLE';
anat = {'HIPPO A'};anat_name = 'HIPPO A';
anat = {'HIPPO M'};anat_name = 'HIPPO M';
anat = {'HIPPO P'};anat_name = 'HIPPO P';
anat = {'PHG'};anat_name = 'PHG';
anat = {'INSULA'};anat_name = 'INSULA';
anat = {'ACC','MCC'};anat_name = 'ACC MCC';
anat = {'ACC'};anat_name = 'ACC';
anat = {'MCC'};anat_name = 'MCC';
anat = {'FG'};anat_name = 'FG';
anat = {'FG','OTS','CS'};anat_name = 'FG OTS CS';
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
anat = {'SFS','SFG'};anat_name = 'SFS SFG';
anat = {'OFC'};anat_name = 'OFC';
anat = {'AMY'};anat_name = 'amygdala';
anat = {'SMA','PARL'};anat_name = 'SMA and PARL';
anat = {'SOG','MOG','IOG'}; anat_name= 'occipital';
anat = {'FOP','POP','TOP'}; anat_name= 'operculum';
anat = {'PCG','CG','POF','CF','PCC'}; anat_name='posterior medial';
anat = {'SPL','IPL','IPS'};anat_name='parieltal area';
anat = {'POSTCENTRAL G','PRECENTRAL G'};anat_name='central area';
anat = {'POSTCENTRAL G'};anat_name='central area';
anat = {'STG','STS','MTG','ITS','ITG'};anat_name='temporal area lateral';
anat = {'PCG'};anat_name='PCG';
anat = {'MTG'};anat_name = 'MTG';
anat = {'IFS','IFG',}
anat = {'INSULA','IFS','IFG','SFS','SFG','ACC','MCC'};anat_name = 'ROI';
anat = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','EC',...
    'FOP','POP','TOP','PARL','INSULA','BASAL'};anat_name='all available';
anat = {'WM','OUT','EMPTY','LESION'};;anat_name='exclude';

side = 'none';%'L','R','none'
site_pick = 'all sites per subj';%all sites per subj

site_pick = 'One site  per subj';%One touch per person


anat = {'STS'};anat_name='STS';
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);


%% choose sites from sheet
% Visit each excel table, add a name column, and concatenate them into a cell
[channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');   
% paticular sites
switch anat_name
    case 'INSULA'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,2,5,6],:);
            T3.anat{3} = T3.anat{3}(3);
        else
        end
    case 'IFS IFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([5,10],:);
            T3.anat{1} = T3.anat{1}(5);
        else
        end
    case 'SFS SFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,4,6],:);
            T3.anat{2} = T3.anat{2}(3);
            T3.anat{3} = T3.anat{3}(1);
        else
        end
end
disp(channame);

%load sheet and vertcat sheet
load('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/mat_files/group_diff_encoding_HFB_reall_behav_sheet.mat')
channame_sheet = cell(length(HFB_behav),1);
for i = 1:length(HFB_behav)
    channame_sheet{i} = HFB_behav{i}.channame{1};
end
HFB_behav = HFB_behav(1,ismember(channame_sheet,channame));
HFB_behav_stats = vertcat(HFB_behav{:});


%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = false;% clusterd based permuation
%%
%make a specific selection of conditions and coloring
if strcmp(project_name,'Grad_CPT')
    conditions = {'mtn','city'};column = 'condNames';
    load cdcol.mat
    plot_params.col = [cdcol.carmine;cdcol.ultramarine];
elseif strcmp(project_name,'race_encoding_simple')
conditions = {'asian','black','white'}; column = 'condNames';
%     conditions = {'own_race','other_races'};column = 'condNames9';
%     conditions = {'my_race_ans','not_my_race_ans'};column = 'condNames8';
%     load cdcol.mat
%     plot_params.col = [cdcol.carmine;cdcol.ultramarine];
end

%%
data_wave = HFB_behav_stats.raw_HFB;
fsample = 500;
winSize = floor(fsample*plot_params.sm);
gusWin = gausswin(winSize)/sum(gausswin(winSize));
data_wave_sm = convn(data_wave,gusWin','same');

recall_isCorrect = HFB_behav_stats.recall_isCorrect;
clean_trial = HFB_behav_stats.clean_trial;
test_race = HFB_behav_stats.test_race;

%concatenate  data for conditions

plot_data_correct = cell(1,length(conditions));
plot_data_incorrect = cell(1,length(conditions));
stats_data_correct = cell(1,length(conditions));
stats_data_incorrect = cell(1,length(conditions));

for i = 1:length(conditions)
    plot_data_correct{i} = data_wave_sm(test_race==i&clean_trial==1&recall_isCorrect==1,:);
    plot_data_incorrect{i} = data_wave_sm(test_race==i&clean_trial==1&recall_isCorrect==0,:);
    stats_data_correct{i} = data_wave(test_race==i&clean_trial==1&recall_isCorrect==1,:);
    stats_data_incorrect{i} = data_wave(test_race==i&clean_trial==1&recall_isCorrect==0,:);
end


range = dsearchn(data_time',[0,1]');
stats_correct = cell(1,length(conditions));
stats_incorrect = cell(1,length(conditions));
for i = 1:length(conditions)
    stats_correct{i} = mean(stats_data_correct{i}(:,range(1):range(2)),2);
    stats_incorrect{i} = mean(stats_data_incorrect{i}(:,range(1):range(2)),2);
end

[h,p,ci,stats] = ttest2(stats_correct{1},stats_incorrect{1});
[h,p,ci,stats] = ttest2(stats_correct{2},stats_incorrect{2});
[h,p,ci,stats] = ttest2(stats_correct{3},stats_incorrect{3});

%% plot figure based on aboving data
clear h
load('cdcol.mat')
figureDim = [100 100 .7 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
cond_names = {'correctly recall','incorrectly recall'};
subtitle = {'Asian','Black','White'};
for ci = 1:length(conditions)
    subplot(1,3,ci)
    
    lineprops.col{1} = plot_params.col(ci,:);
    lineprops.style= '-';
    lineprops.width = plot_params.lw;
    lineprops.edgestyle = '-';%'-';
    hold on
    
    mseb(data_time,nanmean(plot_data_correct{ci}),nanstd(plot_data_correct{ci})/sqrt(size(plot_data_correct{ci},1)),lineprops,1);
    lineprops.style= ':';
    lineprops.col{1} = plot_params.col(ci+5,:);
    mseb(data_time,nanmean(plot_data_incorrect{ci}),nanstd(plot_data_incorrect{ci})/sqrt(size(plot_data_incorrect{ci},1)),lineprops,1); 
    h(1)=plot(data_time,nanmean(plot_data_correct{ci}),'LineWidth',plot_params.lw,'Color',plot_params.col(ci,:));
    h(2)=plot(data_time,nanmean(plot_data_incorrect{ci}),':','LineWidth',plot_params.lw,'Color',plot_params.col(ci+5,:));
    
    xlim(plot_params.xlim)
    
    y_lim = [-.5 1.5];
    ylim_stac = y_lim;
    indx_per = dsearchn(data_time',plot_params.clust_per_win');
    data_time_stac = data_time(indx_per(1):indx_per(2));
    stats_data_correct{ci} = stats_data_correct{ci}(:,indx_per(1):indx_per(2));
    stats_data_incorrect{ci} = stats_data_incorrect{ci}(:,indx_per(1):indx_per(2));
    
    
    data_correct   = permute(stats_data_correct{ci},[3 2 1]);
    data_incorrect = permute(stats_data_incorrect{ci},[3 2 1]);
    rng('default')
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_correct, data_incorrect,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of correct higher than incorrect are ' num2str(clust_info.pos_clust_pval)])
    h(3) = plot(data_time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data_time_stac),'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian lower than black are ' num2str(clust_info.neg_clust_pval)])
    h(4) = plot(data_time_stac, cluster_sig.*ylim_stac(2)-0.05, '-^','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data_time_stac),'MarkerSize',14);
    
    
    
    cond_names_stac = cond_names;
    cond_names_stac{1,1} = 'correctly recall';
    cond_names_stac{1,2} = 'incorrectly recall';
    cond_names_stac{1,3} = 'correct higher';
    cond_names_stac{1,4} = 'correct lower';
    
    
    
    
    
    set(gca,'fontsize',plot_params.textsize)
    box off
    
    
    
    plot([0 0],y_lim, 'Color', [0 0 0], 'LineWidth',2)
    plot(xlim,[0 0], 'Color', [.5 .5 .5], 'LineWidth',1)
    
    
    leg = legend(h,cond_names_stac,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
    legend boxoff
    set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')
    
    xlabel('Time(S)')
    ylabel('HFB Z-score')
    
 
    title(subtitle{ci})
    
end

%% violin plot of the 1sec HFB correct vs incorrect

recall_isCorrect = HFB_behav_stats.recall_isCorrect;
clean_trial = HFB_behav_stats.clean_trial;
test_race = HFB_behav_stats.test_race;
HFB = HFB_behav_stats.("HFB 0.00~1.00");

Asian_c = HFB(recall_isCorrect==1&clean_trial==1&test_race==1);
Asian_ic= HFB(recall_isCorrect==0&clean_trial==1&test_race==1);
Black_c = HFB(recall_isCorrect==1&clean_trial==1&test_race==2);
Black_ic= HFB(recall_isCorrect==0&clean_trial==1&test_race==2);
White_c = HFB(recall_isCorrect==1&clean_trial==1&test_race==3);
White_ic= HFB(recall_isCorrect==0&clean_trial==1&test_race==3);
%%
%z-score outlier exclusion
% Asian_c = Asian_c(abs(zscore(Asian_c))<3);
% Asian_ic= Asian_ic(abs(zscore(Asian_ic))<3);
% Black_c = Black_c(abs(zscore(Black_c))<3);
% Black_ic= Black_ic(abs(zscore(Black_ic))<3);
% White_c = White_c(abs(zscore(White_c))<3);
% White_ic= White_ic(abs(zscore(White_ic))<3);
%% asian
data_stats = [Asian_c;Asian_ic];
group_Asian_c = cell(length(Asian_c),1);
group_Asian_c(:) = deal({'correct'});
group_Asian_ic = cell(length(Asian_ic),1);
group_Asian_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_Asian_c,group_Asian_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','HFB Z-score mean in 1 sec','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(1,:);plot_params.col(6,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(Asian_c,Asian_ic)
[h,p,ci,stats] = ttest2(Asian_c,Asian_ic)
% black
data_stats = [Black_c;Black_ic];
group_Black_c = cell(length(Black_c),1);
group_Black_c(:) = deal({'correct'});
group_Black_ic = cell(length(Black_ic),1);
group_Black_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_Black_c,group_Black_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','HFB Z-score mean in 1 sec','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(2,:);plot_params.col(7,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(Black_c,Black_ic)

% white
data_stats = [White_c;White_ic];
group_White_c = cell(length(White_c),1);
group_White_c(:) = deal({'correct'});
group_White_ic = cell(length(White_ic),1);
group_White_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_White_c,group_White_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','HFB Z-score mean in 1 sec','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(3,:);plot_params.col(8,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(White_c,White_ic)
%% violin plot of the 1sec HFB correct vs incorrect

recall_isCorrect = HFB_behav_stats.recall_isCorrect;
clean_trial = HFB_behav_stats.clean_trial;
test_race = HFB_behav_stats.test_race;
RT = HFB_behav_stats.recall_RT;

Asian_c = RT(recall_isCorrect==1&clean_trial==1&test_race==1);
Asian_ic= RT(recall_isCorrect==0&clean_trial==1&test_race==1);
Black_c = RT(recall_isCorrect==1&clean_trial==1&test_race==2);
Black_ic= RT(recall_isCorrect==0&clean_trial==1&test_race==2);
White_c = RT(recall_isCorrect==1&clean_trial==1&test_race==3);
White_ic= RT(recall_isCorrect==0&clean_trial==1&test_race==3);
%%
%z-score outlier exclusion / need some kind of data trim here?
Asian_c = Asian_c(Asian_c<=3);
Asian_c = Asian_c(abs(zscore(Asian_c))<3);
Asian_ic = Asian_ic(Asian_ic<=3);
Asian_ic= Asian_ic(abs(zscore(Asian_ic))<3);
Black_c = Black_c(Black_c<=3);
Black_c = Black_c(abs(zscore(Black_c))<3);
Black_ic = Black_ic(Black_ic<=3);
Black_ic= Black_ic(abs(zscore(Black_ic))<3);
White_c = White_c(White_c<=3);
White_c = White_c(abs(zscore(White_c))<3);
White_ic = White_ic(White_ic<=3);
White_ic= White_ic(abs(zscore(White_ic))<3);

% asian
data_stats = [Asian_c;Asian_ic];
group_Asian_c = cell(length(Asian_c),1);
group_Asian_c(:) = deal({'correct'});
group_Asian_ic = cell(length(Asian_ic),1);
group_Asian_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_Asian_c,group_Asian_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','RT time','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(1,:);plot_params.col(6,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(Asian_c,Asian_ic)

% black

data_stats = [Black_c;Black_ic];
group_Black_c = cell(length(Black_c),1);
group_Black_c(:) = deal({'correct'});
group_Black_ic = cell(length(Black_ic),1);
group_Black_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_Black_c,group_Black_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','HFB Z-score mean in 1 sec','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(2,:);plot_params.col(7,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(Black_c,Black_ic)

% white
data_stats = [White_c;White_ic];
group_White_c = cell(length(White_c),1);
group_White_c(:) = deal({'correct'});
group_White_ic = cell(length(White_ic),1);
group_White_ic(:) = deal({'incorrect'});%add space this data column will come to the first
data_groups = vertcat(group_White_c,group_White_ic);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','HFB Z-score mean in 1 sec','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',[plot_params.col(3,:);plot_params.col(8,:)]);
g(1,1).axe_property('FontSize',32)


figure('Position',[100 100 600 550]);
g.draw();
[h,p,ci,stats] = ttest2(White_c,White_ic)

