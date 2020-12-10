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



anat = {'STS'};anat_name='STS';
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);


%% Visit each excel table, add a name column, and concatenate them into a cell
[channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');    
%% manully choose one site in each patients
if strcmp(anat_name,'INSULA')
    channame = channame([1,2,5,6],:);
    T3.anat{3} = T3.anat{3}(3);
elseif strcmp(  anat_name,'IFS IFG')
    channame = channame([3,10],:);
    T3.anat{1} = T3.anat{1}(3);
else
end
%% display info of behav data
behv = readtable(['/Users/chao/Desktop/Project_in_Stanford/01_RACE/4_working_data/Behavioral_accuracy/results_summary.xlsx']);
behv_indx = ismember(behv.Chao_patient_ID_in_server,T3.Properties.RowNames);
disp(['the mean of accuracy in ' anat{:} ' is ' num2str(mean(behv.Race_CatAcc_SR(behv_indx))) ' and the std is ' num2str(std(behv.Race_CatAcc_SR(behv_indx)))]);
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
ROL_params.range = [0.1 1];% ROL time range
ROL_params.ratio = 0.5; % percentage of ROL NaN
% ROL-wise HFB
HFB.ROL_onsets = [0 .5];
HFB.ROL_peaks = [.15 .15];

% HFB parameters
tw_series = 1:10;%this is where to adjust
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
             
        end
    else
    end
end

%%
%%stats
% prepare the data
HFB_behav_filt = cell(1,length(channame));
HFB_behav_mix = cell(1,length(channame));

for el = 1:length(channame)
    
    %exclude the HFO and spiky trial
    HFB_behav_filt{el} = HFB_behav{el}(HFB_behav{el}.clean_trial,:);%exclude the HFO and spiky trial
    
    %split the table in to duplicate part and sigle part
    recall_ind = HFB_behav_filt{el}.recall_index;
    [nbin, bin] = histc(recall_ind, unique(recall_ind));
    multiple = find(nbin > 1);
    dup_index    = find(ismember(bin, multiple));
    si_index = setdiff(1:length(recall_ind),dup_index );
    
    si_T = HFB_behav_filt{el}(si_index,:);
    dup_T = HFB_behav_filt{el}(dup_index,:);
    
    %find the duplicate pair in dup_T
    dup_pair = [];
    for i = 1:height(dup_T)
        pair = find(dup_T.recall_index == dup_T.recall_index(i))';
        if pair(1)>pair(2)
            pair = [pair(2),pair(1)];
        else
        end
        dup_pair = [dup_pair;pair];
    end
        dup_pair = unique(dup_pair(:,:),'rows');
        
    %merge duplicate pair
    T_loop1 = [1,2,12:18];% column in the table that is cell or str
    T_loop2 = [3:11,19:size(HFB_behav{1},2)];% column in the table that is number
    
%     v1 = [];
%     for i = 1:length(T_loop1)
%         eval(sprintf('v%d = deal(cell(length(dup_pair),1))', T_loop1(i)))
%     end
%     
%     for i = 1:length(T_loop2)
%         eval(sprintf('v%d = deal([])', T_loop2(i)))
%     end
        

    for k = 1:length(dup_pair)% loop through the duplicated row
        for i = 1:size(HFB_behav{1},2)% loop through the column
            if ismember(i,T_loop1)
                eval(sprintf('v%d(%d,:) = dup_T{:,%d}(dup_pair(%d,1),:);',i,k,i,k))
            elseif ismember(i,T_loop2)
                eval(sprintf('v%d(%d,:) = nanmean(dup_T{:,%d}(dup_pair(%d,:),:));',i,k,i,k))
            else
            end
        end
    end
    
    val_names = [];
    for i = 1:size(HFB_behav{1},2)
         val_names = [val_names sprintf('v%d,',i)];
    end
    val_names(end) = [];
    
    eval(['dup_T_merge = table(' val_names ');']);

    dup_T_merge.Properties.VariableNames = dup_T.Properties.VariableNames;
 
    HFB_behav_mix{el} = vertcat(dup_T_merge, si_T);
end
    HFB_behav_stats = vertcat(HFB_behav_mix{:});
    
    
    
    
    
    
%% set the index
encode_correct = (HFB_behav_stats.isCorrect == 1);

race_asian = (HFB_behav_stats.test_race ==1);
race_black = (HFB_behav_stats.test_race ==2);
race_white = (HFB_behav_stats.test_race ==3);

ROL_indx = (~isnan(HFB_behav_stats.ROL_onsets))&(~isnan(HFB_behav_stats.ROL_peaks));


recall_correct = (HFB_behav_stats.recall_isCorrect == 1);
recall_incorrect = (HFB_behav_stats.recall_isCorrect == 0);

%% comparison between correct and incorrect recall with encoding HFB
idx1 = encode_correct & race_asian & recall_correct;
idx2 = encode_correct & race_asian & recall_incorrect;

idx1 = recall_correct;
idx2 = recall_incorrect;

[h,p,ci,stats] = ttest2(HFB_behav_stats{:,30}(idx1,:),HFB_behav_stats{:,30}(idx2,:))

idx1 = encode_correct & race_black & recall_correct;
idx2 = encode_correct & race_black & recall_incorrect;
[h,p,ci,stats] = ttest2(HFB_behav_stats.HFB_zscore(idx1,:),HFB_behav_stats.HFB_zscore(idx2,:))

idx1 = encode_correct & race_white & recall_correct;
idx2 = encode_correct & race_white & recall_incorrect;
[h,p,ci,stats] = ttest2(HFB_behav_stats.HFB_zscore(idx1,:),HFB_behav_stats.HFB_zscore(idx2,:))
%%
%comparison between correct and incorrect recall with encoding HFB based on
%ROL
idx1 = encode_correct & race_asian & recall_correct & ROL_indx;
idx2 = encode_correct & race_asian & recall_incorrect & ROL_indx;
[h,p,ci,stats] = ttest2(HFB_behav_stats{:,36}(idx1,:),HFB_behav_stats{:,36}(idx2,:))

%%
%correlate the RT vs HFB
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('RT')
end
%correlate the recall_RT vs HFB
idxasian = race_asian;
idxblack = race_black;
idxwhite = race_white;

idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

idxasian = encode_correct & race_asian&ROL_indx;
idxblack = encode_correct & race_black&ROL_indx;
idxwhite = encode_correct & race_white&ROL_indx;


data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats{:,30}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.recall_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1},'type','p');

data_corr_1{1,2} = HFB_behav_stats{:,30}(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.recall_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2},'type','p');

data_corr_1{1,3} = HFB_behav_stats{:,30}(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.recall_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,3},data_corr_2{1,3},'type','p');

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('recall_RT')
end


%correlate the certainty vs HFB
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats{:,31}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.certainty(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1},'type','p');

data_corr_1{1,2} = HFB_behav_stats{:,31}(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.certainty(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2},'type','p');

data_corr_1{1,3} = HFB_behav_stats{:,31}(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.certainty(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2},'type','p');

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('certainty')
end
%correlate the certainty_RT vs HFB
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats{:,31}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.certainty_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats{:,31}(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.certainty_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats{:,31}(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.certainty_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,3},data_corr_2{1,3});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('certainty_RT')
end
%correlate the salience vs HFB
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats{:,31}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.salience(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats{:,31}(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.salience(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats{:,31}(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.salience(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,3},data_corr_2{1,3});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('salience')
end
%correlate the salience_RT vs HFB
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats{:,31}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.salience_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats{:,31}(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.salience_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats{:,31}(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.salience_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,3},data_corr_2{1,3});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('salience_RT')
end
%correlate the HFB vs ROL
idxasian = encode_correct & race_asian&ROL_indx
idxblack = encode_correct & race_black&ROL_indx;
idxwhite = encode_correct & race_white&ROL_indx;

idxasian = race_asian&ROL_indx;
idxblack = race_black&ROL_indx;
idxwhite = race_white&ROL_indx;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats.ROL_onsets(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats{:,36}(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.ROL_onsets(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats{:,36}(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.ROL_onsets(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats{:,36}(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,3},data_corr_2{1,3});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('ROL_HFB')
end

%%
test = vertcat(HFB_behav{:});
race_asian = (test.test_race ==1)&test.clean_trial;
race_black = (test.test_race ==2)&test.clean_trial;
race_white = (test.test_race ==3)&test.clean_trial;
test.ROL_onsets(race_asian,:)
test.ROL_onsets(race_black,:)
test.ROL_onsets(race_white,:)
nanmean(test.ROL_onsets(race_asian,:))
nanmean(test.ROL_onsets(race_black,:))
nanmean(test.ROL_onsets(race_white,:))


data_asian = test.ROL_onsets(race_asian,:);
data_black = test.ROL_onsets(race_black,:);
data_white = test.ROL_onsets(race_white,:);
data_anova = [data_asian;data_black;data_white];
group_asian = repmat({'asian'},size(data_asian,1),1);
group_black = repmat({'black'},size(data_black,1),1);
group_white = repmat({'white'},size(data_white,1),1);
group = [group_asian;group_black;group_white];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_asian),nanstd(data_black ),nanstd(data_white)];
m1 = multcompare(stats1,'ctype','lsd')


%%

encode_correct = (HFB_behav_stats.isCorrect == 1);

race_asian = (HFB_behav_stats.test_race ==1);
race_black = (HFB_behav_stats.test_race ==2);
race_white = (HFB_behav_stats.test_race ==3);

idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;


data_asian = [HFB_behav_stats.recall_RT(idxasian,:) HFB_behav_stats{:,32:end}(idxasian,:)];
data_black = [HFB_behav_stats.recall_RT(idxblack,:) HFB_behav_stats{:,32:end}(idxblack,:)];
data_white = [HFB_behav_stats.recall_RT(idxwhite,:) HFB_behav_stats{:,32:end}(idxwhite,:)];




[dataCovMat,datap] = corr(data_asian);
figure(3), clf
imagesc(dataCovMat)
axis image
set(gca,'clim',[-1 1])
title('Data covariance matrix')
xlabel('xx')
ylabel('??')

yticklabels({'1';'2'})
yticklabels(num2str(HFB.range))

yticklabels({HFB_behav_stats.Properties.VariableNames{22};HFB_behav_stats.Properties.VariableNames{32}})

for i=1:size(datap,1)
    for j=1:size(datap,2)
        text(1,j,num2str(round(datap(1,j),2)),'HorizontalAlignment','center','fontsize',15)
    end
end

[dataCovMat,datap]= corr(data_black);
figure(4), clf
imagesc(dataCovMat)
axis image
set(gca,'clim',[-1 1])
title('Data covariance matrix')
xlabel('??')
ylabel('??')
    for j=1:size(datap,2)
        text(1,j,num2str(round(datap(1,j),2)),'HorizontalAlignment','center','fontsize',15)
    end

[dataCovMat,datap]= corr(data_white);
figure(5), clf
imagesc(dataCovMat)
axis image
set(gca,'clim',[-1 1])
title('Data covariance matrix')
xlabel('??')
ylabel('??')
    for j=1:size(datap,2)
        text(1,j,num2str(round(datap(1,j),2)),'HorizontalAlignment','center','fontsize',15)
    end



data_corr_1{1,1} = HFB_behav_stats{:,31}(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.salience_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

%%
% a clear MATLAB workspace is a clear mental workspace
close all; clear; clc

%% simulate data

% simulation parameters
N = 1000;     % time points
M =   20;     % channels

% time vector (radian units)
t = linspace(0,6*pi,N);

% relationship across channels (imposing covariance)
chanrel = sin(linspace(0,2*pi,M))';

% create the data
data = bsxfun(@times,repmat( sin(t),M,1 ),chanrel);
data = data + randn(M,N);

% two ways of visualizing the multichannel data
figure(1), clf
plot(t,bsxfun(@plus,data,(1:M)'*4))
set(gca,'ytick',[],'xlim',t([1 end]))
xlabel('Time (a.u.)')
ylabel('Channel')

figure(2), clf
imagesc(t,[],data)
xlabel('Time (a.u.)')
ylabel('Channel')
set(gca,'clim',[-1 1]*2)

%% now compute the covariance matrix

% note the size of the output!
dataCovMat = cov(data');

figure(3), clf
imagesc(dataCovMat)
axis image
set(gca,'clim',[-1 1]*.5)
title('Data covariance matrix')
xlabel('??')
ylabel('??')

%% and now the correlation matrix

% note the size of the output!
dataCorrMat = corrcoef(data');

figure(4), clf
imagesc(dataCorrMat)
axis image
set(gca,'clim',[-1 1]*.5)
title('Data correlation matrix')
xlabel('??')
ylabel('??')

% a clear MATLAB workspace is a clear mental workspace
close all; clear; clc

%% the example from the video

% raw correlations
rmg = .7;
rsg = .8;
rms = .9;

% partial correlations
rho_mg_s = (rmg - rsg*rms) / ( sqrt(1-rsg^2)*sqrt(1-rms^2) )
rho_sg_m = (rsg - rmg*rms) / ( sqrt(1-rmg^2)*sqrt(1-rms^2) )

%% now for datasets

N = 76;

% correlated datasets
x1 = linspace(1,10,N)' + randn(N,1);
x2 = x1 + randn(N,1);
x3 = 0*x1 + randn(N,1);


% compute the "raw" correlation matrix
[cormatR,cormatP] = corr([x1 x2 x3]);

% inspect some partial correlations
partialcorr(x2,x3,x1)

%% visualize

figure(1), clf

% raw correlations
subplot(121)
imagesc(cormatR), axis square
set(gca,'clim',[-1 1],'xtick',1:3,'ytick',1:3)
xlabel('Channels'), ylabel('Channels')
title('Raw correlation matrix')

% add text correlation values
for i=1:3
    for j=1:3
        text(i,j,num2str(round(cormatR(i,j),2)),'HorizontalAlignment','center','fontsize',15)
    end
end



% partial correlations
parcormat = partialcorr([x1 x2 x3]);

subplot(122)
imagesc(parcormat), axis square
set(gca,'clim',[-1 1],'xtick',1:3,'ytick',1:3)
xlabel('Channels'), ylabel('Channels')
title('Partial correlation matrix')

% add text correlation values
[I,J] = meshgrid(1:3);
text(I(:),J(:),num2str(round(parcormat(:),2)),'HorizontalAlignment','center','FontSize',15)

%% done.

