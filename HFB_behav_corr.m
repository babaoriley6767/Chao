%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/tony/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/tony/Documents/Stanford/code/lbcn_preproc-master/'))
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
anat_displ =  importdata('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
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
behv = readtable(['/Users/tony/Desktop/Project_in_Stanford/01_RACE/4_working_data/Behavioral_accuracy/results_summary.xlsx']);
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
% HFB parameters
% fixed HFB column
HFB.fix_range = [0 1];                
% buffer HFB column
HFB.buffer_range = [0 .5;.1 .6;.2 .7;.3 .8;.4 .9;.5 1;0.6 0.8;0.7 0.9;0.8 1];   % defalt      [0 .5;.1 .6;.2 .7;.3 .8;.4 .9;.5 1]     
% ROL-wise HFB
HFB.ROL_onsets = [0 .5];
HFB.ROL_peaks = [.15 .15];
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
                % step 1 horizontally concat the recall trialinfo
                HFB_behav{elec} = trialinfo_encode.trialinfo;
                HFB_behav{elec}.encoding_seq = trialinfo_recall_raw.seq.active.race;
                HFB_behav{elec}.recall_index = recall_index;
                HFB_behav{elec}.recall_isCorrect = trialinfo_recall.trialinfo.isCorrect(recall_index,:);
                HFB_behav{elec}.recall_RT = trialinfo_recall.trialinfo.RT(recall_index,:);
                HFB_behav{elec}.certainty = trialinfo_recall.trialinfo.certainty(recall_index,:);
                HFB_behav{elec}.certainty_RT = trialinfo_recall.trialinfo.pa_RT_certainty(recall_index,:);
                HFB_behav{elec}.salience = trialinfo_recall.trialinfo.salience(recall_index,:);
                HFB_behav{elec}.salience_RT = trialinfo_recall.trialinfo.pa_RT_salience(recall_index,:);
                HFB_behav{elec}.clean_trial = clean_trial;
                
                % step 2 horizontally concat the ROL data
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
                
                
                % step 3 horizontally concat the HFB
                % HFB in a fixed range
                HFB_timewindow = HFB.fix_range;
                HFB_range = dsearchn(data_all.time',HFB_timewindow');
                HFB_behav{elec}.HFB = mean(data_all.wave(:,HFB_range(1):HFB_range(2)),2);
                HFB_behav{elec}.Properties.VariableNames{'HFB'} = 'fix_range';
                % HFB from a buffer series
                HFB_range = HFB.buffer_range;
                for h = 1:size(HFB_range,1)
                    HFB_timewindow = HFB.buffer_range(h,:);
                    HFB_range = dsearchn(data_all.time',HFB_timewindow');
                    HFB_behav{elec}.HFB = mean(data_all.wave(:,HFB_range(1):HFB_range(2)),2);
                    HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['buffer range ' num2str(HFB.buffer_range(h,:))];
                end
                 % HFB based on ROL_onsets
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
                HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['ROL onset HFB ' num2str(HFB.ROL_onsets(2))];
                % HFB based on ROL_peaks
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
                HFB_behav{elec}.Properties.VariableNames{'HFB'} = ['ROL peaks HFB from ' num2str(HFB.ROL_onsets(1)) ' to ' num2str(HFB.ROL_onsets(2))];
                
             % step 4 horizontally concat the ID info
                HFB_behav{elec}.channame(1:108,:) = {deal(channame{elec})};
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
    [v01,v02,v03,v04,v05,v06,v07,v08,v09,v10...
     ,v18,v19,v20,v21,v22,v23,v24...
     ,v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,v35,v36,v37,v38,v39,v40]= deal([]);
 
    [v11,v12,v13,v14,v15,v16,v17,v41] = deal(cell(length(dup_pair),1));
    
    for k = 1:length(dup_pair)
        v01(k,:) = dup_T{:,1}(dup_pair(k,1),:);
        v02(k,:) = nanmean(dup_T{:,2}(dup_pair(k,:),:));
        v03(k,:) = nanmean(dup_T{:,3}(dup_pair(k,:),:));
        v04(k,:) = nanmean(dup_T{:,4}(dup_pair(k,:),:));
        v05(k,:) = nanmean(dup_T{:,5}(dup_pair(k,:),:));
        v06(k,:) = nanmean(dup_T{:,6}(dup_pair(k,:),:));
        v07(k,:) = nanmean(dup_T{:,7}(dup_pair(k,:),:));
        v08(k,:) = nanmean(dup_T{:,8}(dup_pair(k,:),:));
        v09(k,:) = nanmean(dup_T{:,9}(dup_pair(k,:),:));
        v10(k,:) = nanmean(dup_T{:,10}(dup_pair(k,:),:));
        v11(k,:) = dup_T{:,11}(dup_pair(k,1),:);
        v12(k,:) = dup_T{:,12}(dup_pair(k,1),:);
        v13(k,:) = dup_T{:,13}(dup_pair(k,1),:);
        v14(k,:) = dup_T{:,14}(dup_pair(k,1),:);
        v15(k,:) = dup_T{:,15}(dup_pair(k,1),:);
        v16(k,:) = dup_T{:,16}(dup_pair(k,1),:);
        v17(k,:) = dup_T{:,17}(dup_pair(k,1),:);
        v18(k,:) = dup_T{:,18}(dup_pair(k,1),:);
        v19(k,:) = nanmean(dup_T{:,19}(dup_pair(k,:),:));
        v20(k,:) = nanmean(dup_T{:,20}(dup_pair(k,:),:));
        v21(k,:) = nanmean(dup_T{:,21}(dup_pair(k,:),:));
        v22(k,:) = nanmean(dup_T{:,22}(dup_pair(k,:),:));
        v23(k,:) = nanmean(dup_T{:,23}(dup_pair(k,:),:));
        v24(k,:) = nanmean(dup_T{:,24}(dup_pair(k,:),:));
        v25(k,:) = nanmean(dup_T{:,25}(dup_pair(k,:),:));
        v26(k,:) = nanmean(dup_T{:,26}(dup_pair(k,:),:));
        v27(k,:) = nanmean(dup_T{:,27}(dup_pair(k,:),:));
        v28(k,:) = nanmean(dup_T{:,28}(dup_pair(k,:),:));
        v29(k,:) = nanmean(dup_T{:,29}(dup_pair(k,:),:));
        v30(k,:) = nanmean(dup_T{:,30}(dup_pair(k,:),:));
        v31(k,:) = nanmean(dup_T{:,31}(dup_pair(k,:),:));
        v32(k,:) = nanmean(dup_T{:,32}(dup_pair(k,:),:));
        v33(k,:) = nanmean(dup_T{:,33}(dup_pair(k,:),:));
        v34(k,:) = nanmean(dup_T{:,34}(dup_pair(k,:),:));
        v35(k,:) = nanmean(dup_T{:,35}(dup_pair(k,:),:));
        v36(k,:) = nanmean(dup_T{:,36}(dup_pair(k,:),:));
        v37(k,:) = nanmean(dup_T{:,37}(dup_pair(k,:),:));
        v38(k,:) = nanmean(dup_T{:,38}(dup_pair(k,:),:));
        v39(k,:) = nanmean(dup_T{:,39}(dup_pair(k,:),:));
        v40(k,:) = nanmean(dup_T{:,40}(dup_pair(k,:),:));
        v41(k,:) = dup_T{:,41}(dup_pair(k,1),:);
    end
    
    dup_T_merge = table(v01,v02,v03,v04,v05,v06,v07,v08,v09,v10,v11,v12...
     ,v13,v14,v15,v16,v17,v18,v19,v20,v21,v22,v23,v24...
     ,v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,v35,v36,v37,v38,v39,v40,v41);

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
h1 = cell(length(HFB.buffer_range),1);
h2 = cell(length(HFB.buffer_range),1);
h3 = cell(length(HFB.buffer_range),1);
p1 = cell(length(HFB.buffer_range),1);
p2 = cell(length(HFB.buffer_range),1);
p3 = cell(length(HFB.buffer_range),1);
for i = 30:38
    [~,j] =find(30:38==i); 

idx1 = encode_correct & race_asian & recall_correct;
idx2 = encode_correct & race_asian & recall_incorrect;
[h1{j},p1{j},ci,stats] = ttest2(HFB_behav_stats{:,i}(idx1,:),HFB_behav_stats{:,i}(idx2,:));

idx1 = encode_correct & race_black & recall_correct;
idx2 = encode_correct & race_black & recall_incorrect;
[h2{j},p2{j},ci,stats] = ttest2(HFB_behav_stats{:,i}(idx1,:),HFB_behav_stats{:,i}(idx2,:));

idx1 = encode_correct & race_white & recall_correct;
idx2 = encode_correct & race_white & recall_incorrect;
[ {j},p3{j},ci,stats] = ttest2(HFB_behav_stats{:,i}(idx1,:),HFB_behav_stats{:,i}(idx2,:));
end
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

data_corr_1 = cell(length(30:35),3);
data_corr_2 = cell(length(30:35),3);
rho = cell(length(30:35),3);
pval = cell(length(30:35),3);
for i = 1: length(30:35) % 2022 JNS
data_corr_1{i,1} = HFB_behav_stats{:,i+29}(idxasian,:);
data_corr_2{i,1} = HFB_behav_stats.recall_RT(idxasian,:);
[rho{i,1},pval{i,1}] = corr(data_corr_1{i,1},data_corr_2{i,1},'type','p');

data_corr_1{i,2} = HFB_behav_stats{:,i+29}(idxblack,:);
data_corr_2{i,2} = HFB_behav_stats.recall_RT(idxblack,:);
[rho{i,2},pval{i,2}] = corr(data_corr_1{i,2},data_corr_2{i,2},'type','p');

data_corr_1{i,3} = HFB_behav_stats{:,i+29}(idxwhite,:);
data_corr_2{i,3} = HFB_behav_stats.recall_RT(idxwhite,:);
[rho{i,3},pval{i,3}] = corr(data_corr_1{i,3},data_corr_2{i,3},'type','p');
end

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

data_corr_1 = cell(length(30:35),3);
data_corr_2 = cell(length(30:35),3);
rho = cell(length(30:35),3);
pval = cell(length(30:35),3);
for i = 1: length(30:35) % 2022 JNS
data_corr_1{i,1} = HFB_behav_stats{:,i+29}(idxasian,:);
data_corr_2{i,1} = HFB_behav_stats.certainty(idxasian,:);
[rho{i,1},pval{i,1}] = corr(data_corr_1{i,1},data_corr_2{i,1},'type','p');

data_corr_1{i,2} = HFB_behav_stats{:,i+29}(idxblack,:);
data_corr_2{i,2} = HFB_behav_stats.certainty(idxblack,:);
[rho{i,2},pval{i,2}] = corr(data_corr_1{i,2},data_corr_2{i,2},'type','p');

data_corr_1{i,3} = HFB_behav_stats{:,i+29}(idxwhite,:);
data_corr_2{i,3} = HFB_behav_stats.certainty(idxwhite,:);
[rho{i,3},pval{i,3}] = corr(data_corr_1{i,3},data_corr_2{i,3},'type','p');
end

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

%correlate the HFB vs accuracy
idxasian = encode_correct & race_asian&ROL_indx;
idxblack = encode_correct & race_black&ROL_indx;
idxwhite = encode_correct & race_white&ROL_indx;

idxasian = race_asian&ROL_indx;
idxblack = race_black&ROL_indx;
idxwhite = race_white&ROL_indx;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats.recall_isCorrect(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats{:,31}(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.recall_isCorrect(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats{:,31}(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.recall_isCorrect(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats{:,31}(idxwhite,:);
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