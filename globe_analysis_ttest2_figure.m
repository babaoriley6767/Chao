%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

% define the cohort
sbj_names_all = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
    ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
    ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
    ;'S17_116_AA';'S17_118_TW';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};

sbj_names = cell(3,1);
sbj_names{1} = sbj_names_all(1:36);%China
sbj_names{2} = sbj_names_all(37:end);%Stanford
sbj_names{3} = sbj_names_all;%all

%make a specific selection of anatomical structures
anat = cell(7,1);
anat_name = cell(7,1);
anat = {'INSULA'};anat_name = 'INSULA';
anat{2} = {'ACC','MCC'};anat_name{2} = 'ACC MCC';
anat{3} = {'FG','OTS','CS'};anat_name{3} = 'FG OTS CS';
anat{4} = {'IFS','IFG'};anat_name{4} = 'IFS IFG';%
anat{5} = {'SFS','SFG'};anat_name{5} = 'SFS SFG';
anat{6} = {'OFC'};anat_name{6} = 'OFC';
anat{7} = {'AMY'};anat_name{7} = 'amygdala';

for iant = 1:7
    for isbj = 1:3
%% Visit each excel table, add a name column, and concatenate them into a cell
T = cell(size(sbj_names,1), 1);
for i = 1:length(sbj_names)
    %cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])%plz adjust accordingly to your ecosystem
    T{i} = readtable(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i} '/' sbj_names{i} '_stats.xlsx']);
    sbj_name_channame = cell(size(T{i},1),1);
    for j = 1:size(T{i},1)
        sbj_name_channame{j} = [sbj_names{i},'-',T{i}.glv_channame{j}];
    end
    T{i}.sbj_name_channame = sbj_name_channame;
end
%Creat another table with rows of specific cohorts and column of specific anatomical
%structures
sz = [size(sbj_names,1) size(anat,2)];
varTypes = cell(1,size(anat,2));
varTypes(:) = {'cell'};
T2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',anat,'RowNames',sbj_names);
%put the glv_index into each space of the table as a vector
for i = 1:length(sbj_names)
    for j = 1:length(anat)
        idx1 = strcmp(T{i}.label,anat{j});
        idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
        idx = idx1 & idx2;
        T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
    end
end
%Since there may be empty electrodes in the stats sheet, Glv_index may be str. Here, all Glv_index in T2 will be transformed into vetor
for i = 1:length(sbj_names)
    for j = 1:length(anat)
        if iscell(T2{sbj_names{i},anat{j}}{:})
            T2{sbj_names{i},anat{j}}{:} = str2double(T2{sbj_names{i},anat{j}}{:});
        end
    end
end
%Creat a third table that horizontally concatenate all the specific
%anatomical structures in together and get rid of the empty rows in T3
sz = [size(sbj_names,1) 1];
varTypes = cell(1,1);
varTypes(:) = {'cell'};
T3 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',{'anat'},'RowNames',sbj_names);
for i = 1:size(T3,1)
    T3{i,:}{:} = horzcat(T2{i,:}{:});
end
loc=cellfun('isempty', T3{:,'anat'} );% 
T3(loc,:)=[];
%%
%define the plot and stats parameters first
project_name = 'race_encoding_simple';
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 20;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
%%
%make a specific selection of conditions
%conditions = {'asian','black','white'}; column = 'condNames';
conditions = {'own_race','other_races'};column = 'condNames9';
if strcmp(column,'condNames9')
    load cdcol.mat
    plot_params.col = [cdcol.ultramarine;cdcol.carmine];
else
end

%%
%concatenate  data for conditions
plot_data = cell(1,length(conditions));
plot_data_all = cell(1,length(conditions));
stats_data = cell(1,length(conditions));
stats_data_all = cell(1,length(conditions));
for i = 1:length(T3.Properties.RowNames)
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
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_name,filesep,'global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
        for j = 1:length(T3.anat{i})
            data_all = concatBlocks(sbj_name, block_names,dirs,T3.anat{i}(j),'HFB','Band',{'wave'},'stimlock_bl_corr');
            
            %smooth is in the trial level
            winSize = floor(data_all.fsample*plot_params.sm);%this part is smooth1
            gusWin = gausswin(winSize)/sum(gausswin(winSize));%this part is smooth2
            data_all.wave_sm = convn(data_all.wave,gusWin','same');%this part is smooth3
            
            [grouped_trials_all,~] = groupConds(conditions,data_all.trialinfo,column,'none',[],false);
            [grouped_trials,cond_names] = groupConds(conditions,data_all.trialinfo,column,stats_params.noise_method,stats_params.noise_fields_trials,false);
            % this part is to exclude HFB over a fixed threshold
            if plot_params.single_trial_replot
                thr_raw =[];
                for di = 1:size(data_all.wave,1)
                    if ~isempty(find(data_all.wave(di,:)>=plot_params.single_trial_thr));
                        fprintf('You have deleted the data over threshold %d from the data \n',plot_params.single_trial_thr);
                    else
                    end
                end
                [thr_raw,thr_column] = find(data_all.wave >= plot_params.single_trial_thr);
                thr_raw = unique(thr_raw);
            end

                for ci = 1:length(conditions)
                    grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);% 
                    plot_data{ci} = [plot_data{ci};nanmean(data_all.wave_sm(grouped_trials{ci},:),1)];
                    plot_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave_sm(grouped_trials_all{ci},:),1)];
                    stats_data{ci} = [plot_data{ci};nanmean(data_all.wave(grouped_trials{ci},:),1)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave(grouped_trials_all{ci},:),1)];
                end
        end
    else
    end
end
%randoml
% 2 sample test
data_win = min(find(data_all.time>stats_params.task_win(1))):max(find(data_all.time<stats_params.task_win(2)));

stats_data{1} = mean(stats_data{1}(:,data_win),2);
stats_data{2} = mean(stats_data{2}(:,data_win),2);

[H,P,CI,STATS] = ttest2(stats_data{1},stats_data{2}); STATS.H = H; STATS.P = P; STATS.CI = CI;
    end
end