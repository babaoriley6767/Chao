%get ready with the addpath,plz adjust accordingly
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
anat = {'INSULA','ACC','MCC','IFS','IFG','SFS','SFG'};anat_name = 'INSULA';
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


%% Visit each excel table, add a name column, and concatenate them into a cell
[channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');    
%% manully choose one site in each patients
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
%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
plot_params.single_trial = false;
plot_params.clust_per = true;% clusterd based permuation
% stats parameters
stats_params = genStatsParams(project_name);

% project and condition info
project_name = 'race_encoding_simple';
conditions = {'asian','black','white'}; column = 'condNames';

% ROL from Jessica and Omri NC parameters
ROL_params = genROLParams_NC(project_name);% set the parameters
ROL_params.range = [0.1 0.8];% ROL time range
ROL_params.ratio = 0.5; % percentage of ROL NaN
% ROL-wise HFB
HFB.ROL_onsets = [0 .5];
HFB.ROL_peaks = [.15 .15];

% HFB parameters
tw_series = 10;%this is where to adjust
tw_HFB = 0:0.1:1;
HFB.range = [];
for i = 1:length(tw_series)
buffer_mitrx = buffer(tw_HFB,tw_series(i)+1,tw_series(i),'nodelay');
HFB_range_loop = [buffer_mitrx(1,:);buffer_mitrx(end,:)]';
HFB.range = [HFB.range; HFB_range_loop];
end
disp('the following is the HFB range we use' );
disp(num2str(HFB.range));
            

%% Horizontally prepare the data
%creat a cell with N tables,each table has the HFB zscore and behav data
HFB_behav = cell(1,length(channame));
for i = 1:length(T3.Properties.RowNames)%loop across patients
    if ~isempty(T3.anat{i})
        indx = i;
        sbj_name = T3.Properties.RowNames{indx};%set basic pipeline parameters
        if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        block_names = BlockBySubj(sbj_name,project_name);
        dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_name,filesep,'global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])%loda globalVar
        
        trialinfo_encode = load(['/Volumes/CHAO_IRON_M/data/psychData/',sbj_name,filesep,block_names{:},filesep,'trialinfo_', block_names{:} '.mat']);%this is the encode behav .mat file
        
        soda_name = dir(fullfile(globalVar.psych_dir, 'sodata*.mat'));
        trialinfo_recall_raw = load([globalVar.psych_dir '/' soda_name.name]);%this is the raw recall behav .mat file
        
        block_name_recall = BlockBySubj(sbj_name,'race_recall');
        trialinfo_recall = load(['/Volumes/CHAO_IRON_M/data/psychData/',sbj_name,filesep,block_name_recall{:},filesep,'trialinfo_', block_name_recall{:} '.mat']);%%this is the recall behav .mat file
        
        [~,recall_index] = ismember(string(trialinfo_recall_raw.seq.active.race),string(trialinfo_recall_raw.seq.memory.active.race));%% key variable that link ecode and recall info, each patient has a unique recall_index
   
        for j = 1:length(T3.anat{i})%loop across site
            
            data_all = concatBlocks(sbj_name, block_names,dirs,T3.anat{i}(j),'HFB','Band',{'wave'},['stimlock_bl_corr']);%'stimlock_bl_corr'
            
               
            [grouped_trials_all,~] = groupConds(conditions,data_all.trialinfo,column,'none',[],false);
            [grouped_trials,cond_names] = groupConds(conditions,data_all.trialinfo,column,stats_params.noise_method,stats_params.noise_fields_trials,false);
           
            % this part is to exclude HFB over a fixed threshold
            if plot_params.single_trial_replot
                thr_raw =[];
                for di = 1:size(data_all.wave,1)
                    if ~isempty(find(data_all.wave(di,:)>=plot_params.single_trial_thr))
                        fprintf('You have deleted the data over threshold %d from the data \n',plot_params.single_trial_thr);
                    else
                    end
                end
                [thr_raw,thr_column] = find(data_all.wave >= plot_params.single_trial_thr);
                thr_raw = unique(thr_raw);
            end
            
            for cond = 1:length(conditions)
                grouped_trials{cond} = setdiff(grouped_trials{cond},thr_raw);
            end
            
            clean_trial = ismember([1:108]',sort(vertcat(grouped_trials{:})));
            
            
            
            if i == 1 % determine which site the code is working at
                elec = j;
            else
                elec = sum(cellfun(@numel,T3.anat(1:(i-1))))+j;
            end
            
                % step 1 horizontally concat the ID info and move it to the
                % second column
                HFB_behav{elec} = trialinfo_encode.trialinfo;
                HFB_behav{elec}.channame(1:108,:) = {deal(channame{elec})};
                HFB_behav{elec}=[HFB_behav{elec}(:,1) HFB_behav{elec}(:,18) HFB_behav{elec}(:,2:17)];
                
                % step 2 horizontally concat the recall trialinfo
                HFB_behav{elec}.encoding_seq = trialinfo_recall_raw.seq.active.race;
                HFB_behav{elec}.recall_index = recall_index;
                HFB_behav{elec}.recall_isCorrect = trialinfo_recall.trialinfo.isCorrect(recall_index,:);
                HFB_behav{elec}.recall_RT = trialinfo_recall.trialinfo.RT(recall_index,:);
                HFB_behav{elec}.certainty = trialinfo_recall.trialinfo.certainty(recall_index,:);
                HFB_behav{elec}.certainty_RT = trialinfo_recall.trialinfo.pa_RT_certainty(recall_index,:);
                HFB_behav{elec}.salience = trialinfo_recall.trialinfo.salience(recall_index,:);
                HFB_behav{elec}.salience_RT = trialinfo_recall.trialinfo.pa_RT_salience(recall_index,:);
                HFB_behav{elec}.clean_trial = clean_trial;
                
                % step 3 horizontally concat the ROL data
                fs = data_all.fsample; % sampling rate
                stimtime = 0;% time 0
                stiminds = find(data_all.time>stimtime,1); % indx of time 0
                befInd = round(ROL_params.pre_event * fs);% indx numbers before time 0
                aftInd = round(ROL_params.dur * fs);% indx numbers after time 0
                data.time = (befInd:aftInd)/fs;% time series
                
                ROL_onsets = nan(108,1);
                ROL_peaks = nan(108,1);
                for cond = 1:length(conditions)
                    data.wave = data_all.wave(grouped_trials{cond},stiminds+befInd:stiminds+aftInd);
                    ROL_data_cond = ROLbootstrap_NC(data, ROL_params);
                    ROL_onsets(grouped_trials{cond}) = ROL_data_cond.onsets';
                    ROL_peaks(grouped_trials{cond}) = ROL_data_cond.peaks';
                end
                
           
                HFB_behav{elec}.ROL_onsets = ROL_onsets;
                HFB_behav{elec}.ROL_peaks = ROL_peaks;
                % Exclude values earlier than 10ms
                HFB_behav{elec}.ROL_onsets(HFB_behav{elec}.ROL_onsets < ROL_params.range(1)) = nan;
%                 HFB_behav{elec}.ROL_onsets(HFB_behav{elec}.ROL_onsets > ROL_params.range(2)) = nan;
                HFB_behav{elec}.ROL_peaks(HFB_behav{elec}.ROL_peaks < ROL_params.range(1)) = nan;
                % ?Sites for which a ROL value could not be obtained in 50% of the trials or more were discarded from the analysis.
                
                
                for cond = 1:length(conditions)
                    idx_race = setdiff(grouped_trials{cond},thr_raw);
                    if strcmp(conditions{cond},'asian') && sum(isnan(HFB_behav{elec}.ROL_onsets(idx_race))) > round(length(HFB_behav{elec}.ROL_onsets(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL onsets in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    elseif strcmp(conditions{cond},'black') && sum(isnan(HFB_behav{elec}.ROL_onsets(idx_race))) > round(length(HFB_behav{elec}.ROL_onsets(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL onsets in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    elseif strcmp(conditions{cond},'white') && sum(isnan(HFB_behav{elec}.ROL_onsets(idx_race))) > round(length(HFB_behav{elec}.ROL_onsets(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL onsets in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    else
                        disp(['there are enought ROL onsets in condition ' conditions{cond} ' in site ' channame{elec,:} ' more than ' num2str(ROL_params.ratio*100) '%'])
                    end
                end
                
                 for cond = 1:length(conditions)
                    idx_race = setdiff(grouped_trials{cond},thr_raw);
                    if strcmp(conditions{cond},'asian') && sum(isnan(HFB_behav{elec}.ROL_peaks(idx_race))) > round(length(HFB_behav{elec}.ROL_peaks(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL peaks  in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    elseif strcmp(conditions{cond},'black') && sum(isnan(HFB_behav{elec}.ROL_peaks(idx_race))) > round(length(HFB_behav{elec}.ROL_peaks(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL peaks  in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    elseif strcmp(conditions{cond},'white') && sum(isnan(HFB_behav{elec}.ROL_peaks(idx_race))) > round(length(HFB_behav{elec}.ROL_peaks(idx_race)))*ROL_params.ratio
                        warning(['not enough ROL peaks  in condition ' conditions{cond} ' in site ' channame{elec,:} ' less than ' num2str(ROL_params.ratio*100) '%'])
                    else
                        disp(['there are enought ROL peaks  in condition ' conditions{cond} ' in site ' channame{elec,:} ' more than ' num2str(ROL_params.ratio*100) '%'])
                    end
                end
               
                
                % step 4.1 horizontally concat HFB based on
                % ROL_onsets/ROL_peaks
                formatSpec = '%.2f';
                for t = 1:108
                    if isnan(HFB_behav{elec}.ROL_onsets(t))
                        HFB_behav{elec}.HFB(t) = NaN;
                    else
                    end
                    HFB_timewindow = zeros(1,2);
                    HFB_timewindow(2) = HFB_behav{elec}.ROL_onsets(t)+HFB.ROL_onsets(2);
                    if HFB_timewindow(2)>ROL_params.range(2)
                        HFB_timewindow(2) = ROL_params.range(2);
                    else
                    end
                    HFB_range = dsearchn(data_all.time',HFB_timewindow');
                    HFB_behav{elec}.HFB(t) = mean(data_all.wave(t,HFB_range(1):HFB_range(2)),2);
                end
                HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['ROL onset HFB ' num2str(HFB.ROL_onsets(1),formatSpec) '~' num2str(HFB.ROL_onsets(2),formatSpec)];
               
                for t = 1:108
                    if isnan(HFB_behav{elec}.ROL_peaks(t))
                        HFB_behav{elec}.HFB(t) = NaN;
                    else
                    end
                    HFB_timewindow = zeros(1,2);
                    HFB_timewindow(1) = HFB_behav{elec}.ROL_peaks(t) - HFB.ROL_peaks(1);
                    HFB_timewindow(2) = HFB_behav{elec}.ROL_peaks(t) + HFB.ROL_peaks(2);   
                    if HFB_timewindow(1)<ROL_params.range(1)
                        HFB_timewindow(1) = ROL_params.range(1);
                    else
                    end
                    if HFB_timewindow(2) > ROL_params.range(2)
                        HFB_timewindow(2) = ROL_params.range(2);
                    else
                    end
                    HFB_range = dsearchn(data_all.time',HFB_timewindow');
                    HFB_behav{elec}.HFB(t) = mean(data_all.wave(t,HFB_range(1):HFB_range(2)),2);
                end
                HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['ROL peaks HFB from -' num2str(HFB.ROL_peaks(1),formatSpec) ' to ' num2str(HFB.ROL_peaks(2),formatSpec)];
                
                
                % step 4.2 horizontally concat the HFB with different range
                for h = 1:size(HFB.range,1)
                    HFB_timewindow = HFB.range(h,:);
                    HFB_range = dsearchn(data_all.time',HFB_timewindow');
                    HFB_behav{elec}.HFB = mean(data_all.wave(:,HFB_range(1):HFB_range(2)),2);
                    HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['HFB ' num2str(HFB.range(h,1),formatSpec) '~' num2str(HFB.range(h,2),formatSpec)];
                end
                
                % step 5 horizontally concat the raw HFB series
                HFB_behav{elec}.raw_HFB = data_all.wave;
             
        end
    else
    end
end
%%
%time-trial-hfb plot of ROL onset asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;


idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);


idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_onset_s,asian_sort] = sort(data_asian_ROL_onset);
data_asian_s = data_asian(asian_sort,:);
data_black_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_onset_s,black_sort] = sort(data_black_ROL_onset);
data_black_s = data_black(black_sort,:);
data_white_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_onset_s,white_sort] = sort(data_white_ROL_onset);
data_white_s = data_white(white_sort,:);

% plot results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_onset_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_onset_s,size(data_black,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_onset_s,size(data_white,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, White'])

data_stats = [data_asian_ROL_onset;data_black_ROL_onset;data_white_ROL_onset];
data_groups = [ones(length(data_asian_ROL_onset),1);ones(length(data_black_ROL_onset),1)*2;ones(length(data_white_ROL_onset),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')% the anova stats of ROL


%%
asian_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_asian,:);
black_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_black,:);
white_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_white,:);



rho = cell(1,3);
pval = cell(1,3);
[rho{1},pval{1}] = corr(data_asian_ROL_onset,asian_hfb);
[rho{2},pval{2}] = corr(data_black_ROL_onset,black_hfb);
[rho{3},pval{3}] = corr(data_white_ROL_onset,white_hfb);



%%
asian_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_asian,:);
black_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_black,:);
white_hfb = HFB_behav_stats.("HFB 0.00~1.00")(idx_clean&idx_ROL_nan&idx_white,:);


ROL_onset_column = data_stats;
HFB_column = [asian_hfb;black_hfb;white_hfb];
group_column = data_groups;

group_column_cell = cell(length(group_column),1);
group_column_cell(group_column==1,1) = deal({'Asian'});
group_column_cell(group_column==2,1) = deal({'Black'});
group_column_cell(group_column==3,1) = deal({'White'});

clear g
close all
figure('Position', [1500 500 700 700]),clf
g(1,1)=gramm('x',ROL_onset_column,'color',group_column_cell);
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','');
% g(1,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',45); %histogram
g(1,1).stat_bin('normalization','cumcount','geom','line')
g(1,1).axe_property('XTickLabel','','Fontsize',15); % We deactivate tht ticks


%Create a scatter plot
g(2,1)=gramm('x',ROL_onset_column,'y',HFB_column,'color',group_column_cell);
g(2,1).set_names('x','ROL onset time','y','HFB zscore','color','Group');
g(2,1).geom_point(); %Scatter plot
g(2,1).stat_glm();
g(2,1).set_point_options('base_size',6);
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend_pos',[0.83 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(2,1).axe_property('Ygrid','on','Fontsize',15); 
g(2,1).legend_axe_handle.FontSize = 25;
% set(g(2,1).legend_axe_handle.FontSize,30);
% h = legend('x','y','z');
%     set(h,'FontSize',28);


%Create y data histogram on the right
g(3,1)=gramm('x',HFB_column,'color',group_column_cell);
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','');
g(3,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',45); %histogram
% g(3,1).stat_bin('normalization','cumcount','geom','stairs')
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');
g(3,1).axe_property('Fontsize',15); 

g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('1s HFB correlate ROL','Fontsize',15);
g.set_color_options('map',plot_params.col);
% g.set_color_options('map','d3_10');
g.draw();
% set(g(2,1).legend_axe_handle.FontSize,30);

%%

load example_data;
load fisheriris.mat

clear g
figure('Position',[100 100 550 550]);

%Create x data histogram on top
g(1,1)=gramm('x',meas(:,2),'color',species);
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','');
g(1,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',15); %histogram
g(1,1).axe_property('XTickLabel',''); % We deactivate tht ticks

%Create a scatter plot
g(2,1)=gramm('x',meas(:,2),'y',meas(:,1),'color',species);
g(2,1).set_names('x','Sepal Width','y','Sepal Length','color','Species');
g(2,1).geom_point(); %Scatter plot
g(2,1).set_point_options('base_size',6);
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend_pos',[0.83 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(2,1).axe_property('Ygrid','on'); 

%Create y data histogram on the right
g(3,1)=gramm('x',meas(:,1),'color',species);
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','');
g(3,1).stat_bin('geom','stacked_bar','fill','all','nbins',15); %histogram
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');

%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Fisher Iris, custom layout');
g.set_color_options('map','d3_10');
g.draw();











%%
% time-trial-hfb plot of ROL peaks asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;



idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_peaks);

idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_peaks_s,asian_sort] = sort(data_asian_ROL_peaks);
data_asian_s = data_asian(asian_sort,:);

data_black_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_peaks_s,black_sort] = sort(data_black_ROL_peaks);
data_black_s = data_black(black_sort,:);

data_white_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_peaks_s,white_sort] = sort(data_white_ROL_peaks);
data_white_s = data_white(white_sort,:);

% plot results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_peaks_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_peaks_s,size(data_black,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_peaks_s,size(data_white,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, White'])

data_stats = [data_asian_ROL_peaks;data_black_ROL_peaks;data_white_ROL_peaks];
data_groups = [ones(length(data_asian_ROL_peaks),1);ones(length(data_black_ROL_peaks),1)*2;ones(length(data_white_ROL_peaks),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')
%%
% time-trial-hfb plot of ROL mixed onset and peaks, asian black and white
%time-trial-hfb plot of ROL onset asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;


idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);


idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_onset_s,asian_sort] = sort(data_asian_ROL_onset);
data_asian_s = data_asian(asian_sort,:);
asian_onsets_peaks = HFB_behav_stats.ROL_peaks(asian_sort,:);

data_black_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_onset_s,black_sort] = sort(data_black_ROL_onset);
data_black_s = data_black(black_sort,:);
black_onsets_peaks = HFB_behav_stats.ROL_peaks(black_sort,:);

data_white_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_onset_s,white_sort] = sort(data_white_ROL_onset);
data_white_s = data_white(white_sort,:);
white_onsets_peaks = HFB_behav_stats.ROL_peaks(white_sort,:);


figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_onset_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
plot(asian_onsets_peaks,size(data_asian,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_onset_s,size(data_black,1):-1:1,'k--','LineWidth',3)
plot(black_onsets_peaks,size(data_black,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_onset_s,size(data_white,1):-1:1,'k--','LineWidth',3)
plot(white_onsets_peaks,size(data_white,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, White'])

data_stats = [data_asian_ROL_onset;data_black_ROL_onset;data_white_ROL_onset];
data_groups = [ones(length(data_asian_ROL_onset),1);ones(length(data_black_ROL_onset),1)*2;ones(length(data_white_ROL_onset),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')

%%
% time-trial-hfb plot of ROL mixed peaks and onsets, asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;



idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_peaks);

idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_peaks_s,asian_sort] = sort(data_asian_ROL_peaks);
data_asian_s = data_asian(asian_sort,:);
asian_peaks_onsets = HFB_behav_stats.ROL_onsets(asian_sort,:);

data_black_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_peaks_s,black_sort] = sort(data_black_ROL_peaks);
data_black_s = data_black(black_sort,:);
black_peaks_onsets = HFB_behav_stats.ROL_onsets(black_sort,:);

data_white_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_peaks_s,white_sort] = sort(data_white_ROL_peaks);
data_white_s = data_white(white_sort,:);
white_peaks_onsets = HFB_behav_stats.ROL_onsets(white_sort,:);

% plot results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_peaks_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
plot(asian_peaks_onsets,size(data_asian,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_peaks_s,size(data_black,1):-1:1,'k--','LineWidth',3)
plot(black_peaks_onsets,size(data_black,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_peaks_s,size(data_white,1):-1:1,'k--','LineWidth',3)
plot(white_peaks_onsets,size(data_white,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, White'])

data_stats = [data_asian_ROL_peaks;data_black_ROL_peaks;data_white_ROL_peaks];
data_groups = [ones(length(data_asian_ROL_peaks),1);ones(length(data_black_ROL_peaks),1)*2;ones(length(data_white_ROL_peaks),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')


