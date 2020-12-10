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


%make a specific selection of cohort
indxcohort = 1:36;%China
sbj_names = sbj_names_all(indxcohort);%China

%make a specific selection of anatomical structures
anat = {'INSULA'};anat_name = 'INSULA';
anat = {'ACC','MCC'};anat_name = 'ACC MCC';
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
anat = {'SFS','SFG'};anat_name = 'SFS SFG';



%% Visit each excel table, add a name column, and concatenate them into a cell
[channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');
%%
%define the plot and stats parameters first
project_name = 'race_encoding_simple';
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 20;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
conditions = {'asian','black','white'}; column = 'condNames';

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
            data_win = dsearchn(data_all.time',stats_params.task_win');
            data_all.wave = data_all.wave(:,data_win(1):data_win(2));
            
            data_stats = [];
            data_groups = [];
            for ci = 1:length(conditions)
                grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);%
                data_stats = [data_stats; mean(data_all.wave(grouped_trials{ci},:),2)];
                data_groups = [data_groups; repmat(ci,length(grouped_trials{ci}),1)];
            end
            
            % ANOVA stats calculation
            [pval,anovatable,stats] = anova1(data_stats,data_groups,'on');
                      
            if i == 1 % determine which site the code is working at
                elec = j;
            else
                elec = sum(cellfun(@numel,T3.anat(1:(i-1))))+j;
            end
            
            channame{elec,2} = anovatable{2,5};
            disp(['the F value of ' channame{elec,1} ' is ' num2str(anovatable{2,5})])
            disp(['And the p value of ' channame{elec,1} ' is ' num2str(pval)])
            disp('the post hoc comparison is as following ')
            multcompare(stats)
            
        end
    end
end
%%
