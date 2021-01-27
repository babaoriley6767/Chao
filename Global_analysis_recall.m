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
% [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');  
stats = 'any_activation';
load('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/cell_of_44_race_cases_tables.mat');%if there is any change of the excel sheet, 
%then this need to update,go to 'Creat_cell_of_tables.mat'
T = T(indxcohort,1);
channame = [];
%Creat another table with rows of specific cohorts and column of specific anatomical
%structures
sz = [size(sbj_names,1) size(anat,2)];
varTypes = cell(1,size(anat,2));
varTypes(:) = {'cell'};
T2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',anat,'RowNames',sbj_names);
%put the glv_index into each space of the table as a vector
if isempty(side)||strcmp(side,'none')   
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            if strcmp(stats,'group_diff')
                idx2 = T{i}.group_diff;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            elseif strcmp(stats,'any_activation')
                idx2 = T{i}.any_activation;
            else
                idx2 = T{i}.all_trials_activation;
            end
            idx = idx1 & idx2;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            channame = [channame; channame_in_T];
        end
    end
else
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            if strcmp(stats,'group_diff')
                idx2 = T{i}.group_diff;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            elseif strcmp(stats,'any_activation')
                idx2 = T{i}.any_activation;
            else
                idx2 = T{i}.all_trials_activation;
            end%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx3 = strcmp(T{i}.LvsR,side);
            idx = idx1 & idx2 & idx3;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
            channame_in_T = T{i}.sbj_name_channame(idx,:);
            channame = [channame; channame_in_T];
        end
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
%% display info of behav data
behv = readtable(['/Users/chao/Desktop/Project_in_Stanford/01_RACE/4_working_data/Behavioral_accuracy/results_summary.xlsx']);
behv_indx = ismember(behv.Chao_patient_ID_in_server,T3.Properties.RowNames);
disp(['the mean of accuracy in ' anat{:} ' is ' num2str(mean(behv.Race_CatAcc_SR(behv_indx))) ' and the std is ' num2str(std(behv.Race_CatAcc_SR(behv_indx)))]);
%%
%define the plot and stats parameters first
project_name ='race_recall';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = false;% clusterd based permuation
%%
%make a specific selection of conditions and coloring

conditions = {'asian','black','white'}; column = 'condNames';


%%
%concatenate  data for conditions
plot_data = cell(1,length(conditions));
plot_data_all = cell(1,length(conditions));
plot_data_trials = cell(1,length(conditions));
plot_data_trials_all = cell(1,length(conditions));
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
                    if ~isempty(find(data_all.wave(di,:)>=plot_params.single_trial_thr));
                        fprintf('You have deleted the data over threshold %d from the data \n',plot_params.single_trial_thr);
                    else
                    end
                end
                [thr_raw,thr_column] = find(data_all.wave >= plot_params.single_trial_thr);
                thr_raw = unique(thr_raw);
            end

                for ci = 1:length(conditions)
                    grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);% make the grouped_trial and thr_raw in together
                    plot_data{ci} = [plot_data{ci};nanmean(data_all.wave_sm(grouped_trials{ci},:),1)];% we use smoothed data for plotting
                    plot_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave_sm(grouped_trials_all{ci},:),1)];
                    stats_data{ci} = [stats_data{ci};nanmean(data_all.wave(grouped_trials{ci},:),1)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [stats_data_all{ci};nanmean(data_all.wave(grouped_trials_all{ci},:),1)];                   
%                     plot_data{ci}  = [plot_data{ci};data_all.wave_sm(grouped_trials{ci},:)];% averaged across all trials
%                     plot_data_all{ci} = [plot_data_all{ci};data_all.wave_sm(grouped_trials_all{ci},:)];
%                     stats_data{ci} = [stats_data{ci};data_all.wave(grouped_trials{ci},:)]; % this part of data were prepared for further comparison (original non-smoothed data)
%                     stats_data_all{ci} = [stats_data_all{ci};data_all.wave(grouped_trials_all{ci},:)];
                end
        end
    else
    end
end
% randomly pick half of the other_races with fixed rng
% if strcmp(column,'condNames9')
%     rng('default');
%     half_index = randsample(size(plot_data{2},1),round(size(plot_data{2},1)/2));
%     plot_data{2} = plot_data{2}(half_index,:);
% else
% end
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
    h(5) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.1, '-^','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of asian lower than white are ' num2str(clust_info.neg_clust_pval)])
    h(5) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.1, '-^','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_black, data_white,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of black higher than white are ' num2str(clust_info.pos_clust_pval)])
    h(6) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.15, '-o','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of ablack lower than white are ' num2str(clust_info.neg_clust_pval)])
    h(6) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.15, '-o','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    
    
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

legend off
% set(gca,'XLabel','Time(S)');%chao
% xlabel('Time(S)') 

sites_num = sum(cellfun(@numel, T3{:,'anat'} ));
sbj_names_num = size(T3,1);
% title([num2str(sites_num),' sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])
disp([num2str(sites_num),' sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])
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
%%
 
data_all.time;
range_baseline = dsearchn(data_all.time',[-.2;0]);
range_HFB = dsearchn(data_all.time',[0;1]);
rng default

dataA = cell(1,length(conditions));
dataB = cell(1,length(conditions));
for ci = 1:length(conditions)
    dataA = mean(stats_data{ci}(:,range_HFB(1):range_HFB(2)),2);
    dataB = mean(stats_data{ci}(:,range_baseline(1):range_baseline(2)),2);
    p = permutation_paired(dataA,dataB,stats_params.nreps);
    disp(['The p value of ' conditions{ci} ' compare to baseline is ' num2str(p)])
end

%%




        stats_params = genStatsParams(project_name);
        stats_params.paired = true;
        stats_params.all_trials = false;
        
        conds = {'asian'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_asian_bl = intersect(find(sig==1),find(greater==1));
        sig_non_asian_bl = find(sig==0);
        
       
        conds = {'black'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_black_bl = intersect(find(sig==1),find(greater==1));
        sig_non_black_bl = find(sig==0);
        
        conds = {'white'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_white_bl = intersect(find(sig==1),find(greater==1));
        sig_non_white_bl = find(sig==0);
        
        ind_abw = union(union(sig_asian_bl,sig_black_bl),sig_white_bl);%%% 

%% plot through each site and make sure of the if there're anyactivation and groupdiff





% group diff
        if contains(sbj_names{i},'C17')||contains(sbj_names{i},'C18')||contains(sbj_names{i},'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        dirs = InitializeDirs(project_name, sbj_names{i}, comp_root, server_root, code_root);
        block_names = BlockBySubj(sbj_names{i},project_name);
        cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_names{i},filesep,'global_',project_name,'_',sbj_names{i},'_',block_names{1},'.mat'])
        T = readtable([sbj_names{i},'_stats.xlsx']);
        
        stats_params = genStatsParams(project_name);
        stats_params.paired = false;
        stats_params.all_trials = false;
        

        conds = {'asian','black'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_asian_black = intersect(find(sig==1),find(greater==1));
        elec_ab2 = globalVar.channame(sig_asian_black);
        sig_black_asian = intersect(find(sig==1),find(greater==0));
        elec_ba2 = globalVar.channame(sig_black_asian);
        
        conds = {'asian','white'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_asian_white = intersect(find(sig==1),find(greater==1));
        elec_aw2 = globalVar.channame(sig_asian_white);
        sig_white_asian = intersect(find(sig==1),find(greater==0));
        elec_wa2 = globalVar.channame(sig_white_asian);
        
        conds = {'black','white'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_black_white = intersect(find(sig==1),find(greater==1));
        elec_bw2 = globalVar.channame(sig_black_white);
        sig_white_black = intersect(find(sig==1),find(greater==0));
        elec_wb2 = globalVar.channame(sig_white_black);
        
        
        vectors = {sig_asian_black, sig_black_asian, sig_asian_white, sig_white_asian, sig_black_white, sig_white_black};  %// row vector
        ind_abw_group_diff = unique(cat(1, vectors{:}));  
               
        
        if ~iscell(T.glv_index)
            group_diff = ismember(T.glv_index,ind_abw_group_diff);
        else
            ind_abw_group_diff_str = cell(length(ind_abw_group_diff),1);
            for j = 1:length(ind_abw_group_diff)
                ind_abw_group_diff_str{j} = num2str(ind_abw_group_diff(j),'%2d');
            end
              group_diff = ismember(T.glv_index,ind_abw_group_diff_str);
            %any_activation = ismember(T.glv_index,mat2cell(num2str(ind_abw,'%2d'),ones(1,length(ind_abw))));
        end
        T.group_diff = group_diff;
        writetable(T, [sbj_names{i} '_stats.xlsx'] );
        
        % any activation
         if contains(sbj_names{i},'C17')||contains(sbj_names{i},'C18')||contains(sbj_names{i},'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        dirs = InitializeDirs(project_name, sbj_names{i}, comp_root, server_root, code_root);
        block_names = BlockBySubj(sbj_names{i},project_name);
        cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_names{i},filesep,'global_',project_name,'_',sbj_names{i},'_',block_names{1},'.mat'])
        T = readtable([sbj_names{i},'.xlsx']);
        
        stats_params = genStatsParams(project_name);
        stats_params.paired = true;
        stats_params.all_trials = false;
        conds = {'asian'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_asian_bl = intersect(find(sig==1),find(greater==1));
        sig_non_asian_bl = find(sig==0);
        
       
        conds = {'black'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_black_bl = intersect(find(sig==1),find(greater==1));
        sig_non_black_bl = find(sig==0);
        
        conds = {'white'};
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr','condNames',conds,'Band','HFB',stats_params);
        sig_white_bl = intersect(find(sig==1),find(greater==1));
        sig_non_white_bl = find(sig==0);
        
        ind_abw = union(union(sig_asian_bl,sig_black_bl),sig_white_bl);%%% 
        
