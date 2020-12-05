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
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = true;% clusterd based permuation
stats_params.HFB_range = [.1 .5]';
%%
project_name = 'race_encoding_simple';
conditions = {'asian','black','white'}; column = 'condNames';

%%
%creat a cell with N tables,each table has the HFB zscore and behav data
HFB_behav = cell(1,length(channame));
for i = 3:length(T3.Properties.RowNames)%loop across patients
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
            
            %smooth is in the trial level
            winSize = floor(data_all.fsample*plot_params.sm);
            gusWin = gausswin(winSize)/sum(gausswin(winSize));
            data_all.wave_sm = convn(data_all.wave,gusWin','same');
            
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
            grouped_trial_after_threshold = setdiff(vertcat(grouped_trials{:}),thr_raw);
            clean_trial = ismember([1:108]',grouped_trial_after_threshold);%
            
            
            
            if i == 1 % determine which site the code is working at
                elec = j;
            else
                elec = sum(cellfun(@numel,T3.anat(1:(i-1))))+j;
            end
                % horizontally concat the encoding and recall trialinfo
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
                
                HFB_range = dsearchn(data_all.time',stats_params.HFB_range);
                HFB_behav{elec}.HFB_zscore = mean(data_all.wave(:,HFB_range(1):HFB_range(2)),2);
                
                 
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
    
    
    HFB_behav_filt{el} = HFB_behav{el}(HFB_behav{el}.clean_trial,:);%exclude the HFO and spiky trial
    
    recall_ind = HFB_behav_filt{el}.recall_index;%split the table in to duplicate part and sigle part
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
     [center,race,gender,StimulusOnsetTime,RT,test_race,test_gender,test_ans,pa_judg,isCorrect...
        ,encoding_seq,recall_index,recall_isCorrect,recall_RT...
        ,certainty,certainty_RT,salience,salience_RT,clean_trial,HFB_zscore] = deal([]);
    [condNames,condNames2,condNames3,condNames4,condNames5,condNames8,condNames9] = deal(cell(length(dup_pair),1));
   
    
    center = dup_T.center(1:length(dup_pair),:);
    race = dup_T.race(1:length(dup_pair),:);
    gender = dup_T.race(1:length(dup_pair),:);
    for k = 1:length(dup_pair)
        StimulusOnsetTime(k,:) = mean(dup_T.StimulusOnsetTime(dup_pair(k,:),:));
        RT(k,:) = mean(dup_T.RT(dup_pair(k,:),:));
        test_race(k,:) = mean(dup_T.test_race(dup_pair(k,:),:));
        test_gender(k,:) = mean(dup_T.test_gender(dup_pair(k,:),:));
        test_ans(k,:) = mean(dup_T.test_ans(dup_pair(k,:),:));
        pa_judg(k,:) = mean(dup_T.pa_judg(dup_pair(k,:),:));
        isCorrect(k,:) = mean(dup_T.isCorrect(dup_pair(k,:),:));
        condNames(k,:) = dup_T.condNames(dup_pair(k,1),:);
        condNames2(k,:) = dup_T.condNames2(dup_pair(k,1),:);
        condNames3(k,:) = dup_T.condNames3(dup_pair(k,1),:);
        condNames4(k,:) = dup_T.condNames4(dup_pair(k,1),:);
        condNames5(k,:) = dup_T.condNames5(dup_pair(k,1),:);
        condNames8(k,:) = dup_T.condNames8(dup_pair(k,1),:);
        condNames9(k,:) = dup_T.condNames9(dup_pair(k,1),:);
        encoding_seq(k,:) = dup_T.encoding_seq(dup_pair(k,1),:);
        recall_index(k,:) = mean(dup_T.recall_index(dup_pair(k,:),:));
        recall_isCorrect(k,:) = mean(dup_T.recall_isCorrect(dup_pair(k,:),:));
        recall_RT(k,:) = mean(dup_T.recall_RT(dup_pair(k,:),:));
        certainty(k,:) = mean(dup_T.certainty(dup_pair(k,:),:));
        certainty_RT(k,:) = mean(dup_T.certainty_RT(dup_pair(k,:),:));
        salience(k,:) = mean(dup_T.salience(dup_pair(k,:),:));
        salience_RT(k,:) = mean(dup_T.salience_RT(dup_pair(k,:),:));
        clean_trial(k,:) = mean(dup_T.clean_trial(dup_pair(k,:),:));
        HFB_zscore(k,:) = mean(dup_T.HFB_zscore(dup_pair(k,:),:));
    end
    dup_T_merge = table(center,race,gender,StimulusOnsetTime,RT,test_race,test_gender,test_ans,pa_judg,isCorrect...
        ,condNames,condNames2,condNames3,condNames4,condNames5,condNames8,condNames9,encoding_seq,recall_index,recall_isCorrect,recall_RT...
        ,certainty,certainty_RT,salience,salience_RT,clean_trial,HFB_zscore);
    
    HFB_behav_mix{el} = vertcat(dup_T_merge, si_T);
    
    
end


HFB_behav_stats = vertcat(HFB_behav_mix{:});
%1 inclusion
%% set the index
encode_correct = (HFB_behav_stats.isCorrect == 1);

race_asian = (HFB_behav_stats.test_race ==1);
race_black = (HFB_behav_stats.test_race ==2);
race_white = (HFB_behav_stats.test_race ==3);

%% comparison between correct and incorrect recall with encoding HFB
recall_correct = (HFB_behav_stats.recall_isCorrect == 1);
recall_incorrect = (HFB_behav_stats.recall_isCorrect == 0);

idx1 = encode_correct & race_asian & recall_correct;
idx2 = encode_correct & race_asian & recall_incorrect;
[h,p,ci,stats] = ttest2(HFB_behav_stats.HFB_zscore(idx1,:),HFB_behav_stats.HFB_zscore(idx2,:))

idx1 = encode_correct & race_black & recall_correct;
idx2 = encode_correct & race_black & recall_incorrect;
[h,p,ci,stats] = ttest2(HFB_behav_stats.HFB_zscore(idx1,:),HFB_behav_stats.HFB_zscore(idx2,:))

idx1 = encode_correct & race_white & recall_correct;
idx2 = encode_correct & race_white & recall_incorrect;
[h,p,ci,stats] = ttest2(HFB_behav_stats.HFB_zscore(idx1,:),HFB_behav_stats.HFB_zscore(idx2,:))
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
idxasian = encode_correct & race_asian;
idxblack = encode_correct & race_black;
idxwhite = encode_correct & race_white;

data_corr_1 = cell(1,3);
data_corr_2 = cell(1,3);
rho = cell(1,3);
pval = cell(1,3);

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.recall_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.recall_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.recall_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

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

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.certainty(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.certainty(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.certainty(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

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

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.certainty_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.certainty_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.certainty_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

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

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.salience(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.salience(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.salience(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

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

data_corr_1{1,1} = HFB_behav_stats.HFB_zscore(idxasian,:);
data_corr_2{1,1} = HFB_behav_stats.salience_RT(idxasian,:);
[rho{1},pval{1}] = corr(data_corr_1{1,1},data_corr_2{1,1});

data_corr_1{1,2} = HFB_behav_stats.HFB_zscore(idxblack,:);
data_corr_2{1,2} = HFB_behav_stats.salience_RT(idxblack,:);
[rho{2},pval{2}] = corr(data_corr_1{1,2},data_corr_2{1,2});

data_corr_1{1,3} = HFB_behav_stats.HFB_zscore(idxwhite,:);
data_corr_2{1,3} = HFB_behav_stats.salience_RT(idxwhite,:);
[rho{3},pval{3}] = corr(data_corr_1{1,2},data_corr_2{1,2});

figure(1),clf
for fi = 1:3
subplot(1,3,fi)
plot(data_corr_1{1,fi},data_corr_2{1,fi},'kp')
xlabel('HFB'),ylabel('salience_RT')
end


%% plot figure based on aboving data
clear h
load('cdcol.mat')
figureDim = [100 100 .23 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
for ci = 1:length(conditions)
    lineprops.col{1} = plot_params.col(ci,:);
    if ~strcmp(plot_params.eb,'none')
        lineprops.style= '-';
        lineprops.width = plot_params.lw;
        lineprops.edgestyle = '-';
        if strcmp(plot_params.eb,'ste')
            mseb(data_all.time,nanmean(plot_data{ci}),nanstd(plot_data{ci})/sqrt(size(plot_data{ci},1)),lineprops,1);
        
            hold on
        else %'std'
            mseb(data_all.time,nanmean(plot_data{ci}),nanstd(plot_data{ci}),lineprops,1);
            hold on
        end
    end
 
    h(ci)=plot(data_all.time,nanmean(plot_data{ci}),'LineWidth',plot_params.lw,'Color',plot_params.col(ci,:));
    hold on
end
xlim(plot_params.xlim)


if contains(sbj_names{1},'C17_20') && contains(sbj_names{end},'S20_152_HT')
    xlabel(plot_params.xlabel);
else
end


if strcmp(anat_name,'INSULA')
    ylabel(plot_params.ylabel);
else
end

set(gca,'fontsize',plot_params.textsize)
box off


    y_lim = ylim;

    % this part is clusterd based permuation
   
if plot_params.clust_per
    
    ylim_stac = y_lim;
    indx_per = find(abs(data_all.time-plot_params.clust_per_win(1))<0.001):find(abs(data_all.time-plot_params.clust_per_win(2))<0.001);
    data.time_stac = data_all.time(indx_per);
    
    data_asian = stats_data{1}(:,indx_per);
    data_asian = permute(data_asian,[3 2 1]);
    data_black = stats_data{2}(:,indx_per);
    data_black = permute(data_black,[3 2 1]);
    data_white = stats_data{3}(:,indx_per);
    data_white = permute(data_white,[3 2 1]);
    rng('default')
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_asian, data_black,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian higher than black are ' num2str(clust_info.pos_clust_pval)])
    h(4) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian lower than black are ' num2str(clust_info.neg_clust_pval)])
    h(4) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_asian, data_white,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian higher than white are ' num2str(clust_info.pos_clust_pval)])
    h(5) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.1, '-^','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerFaceColor',cdcol.chinese_white,'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian lower than white are ' num2str(clust_info.neg_clust_pval)])
    h(5) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.1, '-^','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerFaceColor',cdcol.chinese_white,'MarkerSize',14);
    
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_black, data_white,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of black higher than white are ' num2str(clust_info.pos_clust_pval)])
    h(6) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.15, '-o','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerFaceColor',cdcol.chinese_white,'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of ablack lower than white are ' num2str(clust_info.neg_clust_pval)])
    h(6) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.15, '-o','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerFaceColor',cdcol.chinese_white,'MarkerSize',14);
    
    
    
    cond_names_stac = cond_names;
    cond_names_stac{1,1} = 'asian';
    cond_names_stac{1,2} = 'black';
    cond_names_stac{1,3} = 'white';
    cond_names_stac{1,4} = 'asian vs black';
    cond_names_stac{1,5} = 'asian vs white';
    cond_names_stac{1,6} = 'black vs white';
else
end


% take care of the axis, xlim, ylim, legend etc.
if size(data_all.trialinfo.allonsets,2) > 1
    time_events = cumsum(nanmean(diff(data_all.trialinfo.allonsets,1,2)));
    for i = 1:length(time_events)
        plot([time_events(i) time_events(i)],y_lim,'Color', [.5 .5 .5], 'LineWidth',1)
    end
else
end
plot([0 0],y_lim, 'Color', [0 0 0], 'LineWidth',2)
plot(xlim,[0 0], 'Color', [.5 .5 .5], 'LineWidth',1)
ylim(y_lim)
box on 
leg = legend(h,cond_names,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')

%legend off
% set(gca,'XLabel','Time(S)');%chao
xlabel('Time(S)') 

sites_num = sum(cellfun(@numel, T3{:,'anat'} ));
sbj_names_num = size(T3,1);
title([num2str(sites_num),' sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])
%% anova and post hoc

data_asian = mean(stats_data{1}(:,indx_per),2);
data_black = mean(stats_data{2}(:,indx_per),2);
data_white = mean(stats_data{3}(:,indx_per),2);
data_anova = [data_asian;data_black;data_white];
group_asian = repmat({'asian'},size(data_asian,1),1);
group_black = repmat({'black'},size(data_black,1),1);
group_white = repmat({'white'},size(data_white,1),1);
group = [group_asian;group_black;group_white];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [std(data_asian),std(data_black ),std(data_white)];
m1 = multcompare(stats1,'ctype','tukey-kramer')
% m1 = multcompare(stats1,'ctype','bonferroni')
% m1 = multcompare(stats1,'ctype','lsd')
%% stats
load('cdcol.mat')
figureDim = [100 100 .23 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
for ci = 1:length(conditions)
       plot_params.single_trial
            subplot(3,1,ci)
            plot(data_all.time,plot_data_trials{ci}', 'Color',plot_params.col(ci,:))
            hold on
            y_lim = [-1,10];
            xlim(plot_params.xlim)
            ylim(y_lim)
      
end
%% stats
[H,P,CI,STATS] = ttest2(stats_data{1},stats_data{2}); STATS.H = H; STATS.P = P; STATS.CI = CI;

%% plot the distribution of sites among cases
figureDim = [100 100 .23 .35 ];
figure('units', 'normalized', 'outerposition', figureDim)
%x = 1:sbj_names_num;
x = cell(1,size(T3,1));
for i = 1:size(T3,1)
    x{i} = T3.Row{i};
end
x=categorical(x);
y = cellfun(@numel, T3{:,'anat'} );
barh(x,y)
set(gca,'fontsize',plot_params.textsize)
sbj_names_num = size(T3,1);
title(['distribution of sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])
set(gca,'XTick',[0:max(y)]);
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
%% plot the sites in MNI space 
%This part can work but the figure is not pretty, I'm still working on this part, if you know some way that we could plot a nicer brain
%plz tell me know, chao
%PlotCoverageElect
%PlotCoverageElectComb

%define the figure parameters first
cfg = [];
% cfg.highlight_col = [1 0 0];
% cfg.chan_highlight = [];
% cfg.correction_factor = 0;
% cfg.MarkerSize = 20;
% cfg.alpha = 0.3;
% cfg.lobe = 'medial';
% cfg.correction_factor = 0;
cfg.flip_right = true;% this is to flipped all the sites

% get the MNI coords and L/R information from variable T
coords = struct;
coords.MNI_coord = [];
coords.LvsR = [];
coords.isleft = [];
coords.channame = [];

if isempty(side)||strcmp(side,'none')
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.group_diff;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation/T{i}.group_diff default: any_activation
            idx = idx1 & idx2;
            coords_in_T = [T{i}.MNI_coord_1 T{i}.MNI_coord_2 T{i}.MNI_coord_3];
%             coords_in_T = [T{i}.fsaverageINF_coord_1 T{i}.fsaverageINF_coord_2 T{i}.fsaverageINF_coord_3];
            coords_in_T = coords_in_T(idx,:);
            LvsR_in_T = T{i}.LvsR(idx,:);
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            coords.MNI_coord = [coords.MNI_coord;coords_in_T];
            coords.LvsR = [coords.LvsR;LvsR_in_T];
            coords.channame = [coords.channame; channame_in_T];
        end
    end
else
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx3 = strcmp(T{i}.LvsR,side);
            idx = idx1 & idx2 & idx3;
            coords_in_T = [T{i}.MNI_coord_1 T{i}.MNI_coord_2 T{i}.MNI_coord_3];
            coords_in_T = coords_in_T(idx,:);
            LvsR_in_T = T{i}.LvsR(idx,:);
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            coords.MNI_coord = [coords.MNI_coord;coords_in_T];
            coords.LvsR = [coords.LvsR;LvsR_in_T];
            coords.channame = [coords.channame; channame_in_T];
        end
    end
end


coords.isleft = [];
for i = 1:length(coords.LvsR)
    if strcmp(coords.LvsR{i},'L')
        coords.isleft(i,1) = 1;
    else
        coords.isleft(i,1) = 0;
    end
end


% this is for flip
if cfg.flip_right
for i = 1:size(coords.MNI_coord,1)
    if coords.MNI_coord(i,1) <= 0
        coords.MNI_coord(i,1) = coords.MNI_coord(i,1)+2*abs(coords.MNI_coord(i,1));
    end
end
else
end
coords.isleft = zeros(length(coords.isleft),1);

% plot the brain and the sites

% PlotCoverageGroup(coords,cfg)
%% plot the sites in MNI space by using iELVis
addpath(genpath('/Users/chao/Desktop/function_tools/for_plot/CCEP_fMRI/'))
addpath(genpath('/Users/chao/Desktop/function_tools/for_plot/iELVis-master/'))
global globalFsDir;
globalFsDir ='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';
fsDir='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';
cd([fsDir]);
%plotting=0;
load('DDS_parc_colors.mat');
load('cdcol_stanford.mat');
cdcol.loc_yellow = hex2rgb('#fff301');


% asian race orange #ff7c01     [1 0.4863 0.0039]
% Black race purple #612696     [0.3804 0.1490 0.5882]
% white race green  #02a755     [0.0078 0.6549 0.3333]
% localization yellow #fff301   [1 0.9529 0.0039]




eleColors =[];

for i = 1:size(coords.channame)
%      if contains(coods.channame{i},'C17')|| contains(coods.channame{i},'C18')||contains(coods.channame{i},'C19')
%         eleColors(i,:) = [0.9216,0.5686,0.3059];
% %     elseif flag(i) == 2
% %         eleColors(i,:) = [0.9196,0.1412,0.3451]
%     else
        eleColors(i,:) =[cdcol.loc_yellow];
end

cfg=[];
cfg.view='r';
cfg.elecSize=12;
cfg.surfType='inflated';    
cfg.opaqueness=1;
cfg.ignoreDepthElec='n';
cfg.elecCoord=[coords.MNI_coord coords.isleft];
cfg.elecNames=coords.channame;
cfg.elecColors = eleColors;
cfg.elecColorScale=[0 1];
cfg.title = [];
cfg.backgroundColor = [1,1,1];
cfgOut=plotPialSurf_v2('fsaverage',cfg);




%% Check whether the left and right coordinates are correct in the the subjvar 
for i = 1:length(sbj_names_all)
    load(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/',sbj_names_all{i},'/','subjVar_',sbj_names_all{i},'.mat'])
    elinfo_link = subjVar.elinfo;
    
    for j = 1:size(elinfo_link,1)
        if strcmp(elinfo_link.LvsR{j},'L')&&elinfo_link.MNI_coord(j,1)<0
            disp('good to go')
        elseif strcmp(elinfo_link.LvsR{j},'R')&&elinfo_link.MNI_coord(j,1)>0
            disp('good to go')
        elseif strcmp(elinfo_link.LvsR{j},'L')&&elinfo_link.MNI_coord(j,1)>0
            warning(['please check the subjvar' sbj_names{i}])
            return
        elseif strcmp(elinfo_link.LvsR{j},'R')&& elinfo_link.MNI_coord(j,1)<0
            warning(['please check the subjvar' sbj_names{i}])
            return
        else
            disp('catched empty site')
        end
    end
end
%%





%% pedro method
task = 'MMR'
cond_names = {'math', 'autobio'};
column = 'condNames';
data = concatenate_multiple_elect(elect_list, task, dirs, 'Band', 'HFB', 'stim');
plot_group_elect(data, task, cond_names, column);
%% about the sheet
% take care of the str and num of FS_ind and glv_index
for i = 1:size(T,1)
    if ~iscell(T{i}.FS_ind)
       T{i}.FS_ind = num2cell(T{i}.FS_ind);
       for j = 1:size(T{i},1)
           T{i}.FS_ind{j} = num2str(T{i}.FS_ind{j});
       end
    else
    end 
end
for i = 1:size(T,1)
    if ~iscell(T{i}.glv_index)
       T{i}.glv_index = num2cell(T{i}.glv_index);
       for j = 1:size(T{i},1)
           T{i}.glv_index{j} = num2str(T{i}.glv_index{j});
       end
    else
    end
end

for i = 1:size(T,1)
    if ~iscell(T{i}.label)
       T{i}.label = num2cell(T{i}.label);
       for j = 1:size(T{i},1)
           T{i}.label{j} = num2str(T{i}.label{j});
       end
    else
    end
end
% concatenate table vertically
bigtable = [T{1};T{2}];
for tidx = 3:numel(T)
   bigtable = [bigtable; T{tidx}];
end

%% Pedro's code 

% % Group analyses
% 
% 1. Select which electrodes to include
% -Statistical
%     Compare all trails agains baseline
%     permutation test between avg baseline period whithin trial vs. avg 1s period within trial
% -Anatomical
%     ROIS
% After these steps you should have the electrode_list
% 
% 2. Concatenate all electrodes of interest
chan_num = [];
sbj_name_p = [];
for i = 1:size(T3,1)
   for j = 1:size(T3.anat{i},2)
    a = T3.anat{i}(j);
    chan_num = [chan_num;a];
    b = T3.Row(i);
    sbj_name_p = [sbj_name_p;b];
   end
end

project_name = 'race_encoding_simple';
electrode_list = table(chan_num,sbj_name_p);
dirs = InitializeDirs(project_name, sbj_name_p{1}, comp_root, server_root, code_root);
data_all = concatenate_multiple_elect(electrode_list, project_name, dirs, 'Band', 'HFB', 'stim');
cond_names = {'asian','black','white'};column = 'condNames';
plot_group_elect(data_all, project_name, cond_names, column);
% 3. Compare two conditions across electrodes
% -Define conditions and time window
%     cond_names = {'autobio', 'math'};
%     column = 'condNames';
stats_params = genStatsParams(project_name);
% -Average each electrode per condition within the 1s period
%     (now you have a single value per electrode per condition)
%     Ready to compare electrodes with independent sample t-test
STATS = stats_group_elect(data_all,data_all.time, project_name,cond_names, column, stats_params);

