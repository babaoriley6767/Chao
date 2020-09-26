
%% this is how we determined the sample of C18_37 for ifs_ifg

data =  cellfun('length',T3.anat)'; %this data from global_analysis

data(7) = 7;% this is the largest cite number that can fit the modified Z-scored method

dataMed = median(data);
dataMAD = mad(data,1); % IMPORTANT: note the second input

dataMz = norminv(.75)*(data-dataMed) ./ dataMAD;

zscorethresh = 3
% show the data!
figure(3), clf
subplot(211), hold on
plot(data,'k^','markerfacecolor','w','markersize',12,'LineWidth',2);
set(gca,'xtick',[]), box off
ylabel('Orig. scale')

subplot(212), hold on
plot(dataMz,'k^','markerfacecolor','w','markersize',12,'LineWidth',2);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r','LineWidth',2)
set(gca,'xtick',[]), box off
ylabel('Median dev. units (Mz)')


%%
%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

sbj_names_all = {'C18_37'};
sbj_names = sbj_names_all;
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
side = 'none';%'L','R','none'





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
if isempty(side)||strcmp(side,'none')
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx = idx1 & idx2;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
        end
    end
else
    for i = 1:length(sbj_names)
        for j = 1:length(anat)
            idx1 = strcmp(T{i}.label,anat{j});
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
            idx3 = strcmp(T{i}.LvsR,side);
            idx = idx1 & idx2 & idx3;
            T2{sbj_names{i},anat{j}} = {T{i}.glv_index(idx)'};
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
%%
%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 20;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = true;
%%
%make a specific selection of conditions and coloring
if strcmp(project_name,'Grad_CPT')
    conditions = {'mtn','city'};column = 'condNames';
    load cdcol.mat
    plot_params.col = [cdcol.carmine;cdcol.ultramarine];
elseif strcmp(project_name,'race_encoding_simple')
    conditions = {'asian','black','white'}; column = 'condNames';
%     conditions = {'own_race','other_races'};column = 'condNames9';
%     load cdcol.mat
%     plot_params.col = [cdcol.carmine;cdcol.ultramarine];
end

%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 20;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = true;
%%
%make a specific selection of conditions and coloring
if strcmp(project_name,'Grad_CPT')
    conditions = {'mtn','city'};column = 'condNames';
    load cdcol.mat
    plot_params.col = [cdcol.carmine;cdcol.ultramarine];
elseif strcmp(project_name,'race_encoding_simple')
    conditions = {'asian','black','white'}; column = 'condNames';
%     conditions = {'own_race','other_races'};column = 'condNames9';
%     load cdcol.mat
%     plot_params.col = [cdcol.carmine;cdcol.ultramarine];
end


%% randomize
tic
n = 10000;%the number of loop
survive = 0; % the number of survive

while n >1
    n = n-1;
T4 = T3;
% s = rng;
rand37 = randperm(24,7);
T4.anat{7} = T4.anat{7}(rand37);

%%
%concatenate  data for conditions
plot_data = cell(1,length(conditions));
plot_data_all = cell(1,length(conditions));
plot_data_trials = cell(1,length(conditions));
plot_data_trials_all = cell(1,length(conditions));
stats_data = cell(1,length(conditions));
stats_data_all = cell(1,length(conditions));
for i = 1:length(T4.Properties.RowNames)
    if ~isempty(T4.anat{i})
        indx = i;
        sbj_name = T4.Properties.RowNames{indx};%set basic pipeline parameters
        if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        block_names = BlockBySubj(sbj_name,project_name);
        dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_name,filesep,'global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
        for j = 1:length(T4.anat{i})
            data_all = concatBlocks(sbj_name, block_names,dirs,T4.anat{i}(j),'HFB','Band',{'wave'},['stimlock_bl_corr']);%'stimlock_bl_corr'
            
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
                    grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);% make the grouped_trial and thr_raw in together
                    stats_data{ci} = [stats_data{ci};nanmean(data_all.wave_sm(grouped_trials{ci},:),1)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [stats_data_all{ci};nanmean(data_all.wave_sm(grouped_trials_all{ci},:),1)];
                end
        end
    else
    end
end
%%

   
    indx_per = find(abs(data_all.time-plot_params.clust_per_win(1))<0.001):find(abs(data_all.time-plot_params.clust_per_win(2))<0.001);
   
    
    data_asian = stats_data{1}(:,indx_per);
    data_asian = permute(data_asian,[3 2 1]);
    data_black = stats_data{2}(:,indx_per);
    data_black = permute(data_black,[3 2 1]);
    data_white = stats_data{3}(:,indx_per);
    data_white = permute(data_white,[3 2 1]);
    
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_asian, data_black,[]);
    cluster_indx_ab = find(clust_info.pos_clust_pval<0.05);
    cluster_sig_ab = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig_ab(cluster_sig==0)=NaN;
    
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_asian, data_white,[]);
    cluster_indx_aw = find(clust_info.pos_clust_pval<0.05);
    cluster_sig_aw = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig_aw(cluster_sig==0)=NaN;
    
    if ~isempty(cluster_indx_ab)||~isempty(cluster_indx_aw)
        survive = survive+1;
    else
    end

end

disp(['the survive ratio of PT37 is ' num2str((survive/10000)*100) ' %']);
toc


%tic
%b = nchoosek(5,2);
%c = nchoosek([1 2 3 4 5],2);
%toc