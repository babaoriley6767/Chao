%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

% define the cohort
sbj_names_all = {'S20_158'};

%make a specific selection of cohort
indxcohort = 1;%China

sbj_names = sbj_names_all(indxcohort);%China


channame = {'S20_158-7L3'};
anat = {[81]};
T3 = table(anat);
T3.Properties.RowNames = {'S20_158'};

%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';%
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = true;% clusterd based permuation
%%
%make a specific selection of conditions and coloring
conditions = {'asian','black','white'}; column = 'condNames';
load cdcol.mat
cdcol.own_race_red =           [0.8200 0      0.1800];  %#D1002E
cdcol.other_race_black =       [0      0.1608 0.2353];  %#00293C
plot_params.col = [cdcol.own_race_red;cdcol.other_race_black];

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
        block_names = block_names(1,3);%%%%%%%%%%%%%%%%% this is where to change
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
%                     plot_data{ci} = [plot_data{ci};nanmean(data_all.wave_sm(grouped_trials{ci},:),1)];% we use smoothed data for plotting
%                     plot_data_all{ci} = [plot_data_all{ci};nanmean(data_all.wave_sm(grouped_trials_all{ci},:),1)];
%                     stats_data{ci} = [stats_data{ci};nanmean(data_all.wave(grouped_trials{ci},:),1)]; % this part of data were prepared for further comparison (original non-smoothed data)
%                     stats_data_all{ci} = [stats_data_all{ci};nanmean(data_all.wave(grouped_trials_all{ci},:),1)];                   
                    plot_data{ci}  = [plot_data{ci};data_all.wave_sm(grouped_trials{ci},:)];% averaged across all trials
                    plot_data_all{ci} = [plot_data_all{ci};data_all.wave_sm(grouped_trials_all{ci},:)];
                    stats_data{ci} = [stats_data{ci};data_all.wave(grouped_trials{ci},:)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [stats_data_all{ci};data_all.wave(grouped_trials_all{ci},:)];
                end
        end
    else
    end
end



plot_data_SRvsOR = cell(1,2);
plot_data_all_SRvsOR = cell(1,2);
stats_data_SRvsOR = cell(1,2);
stats_data_all_SRvsOR = cell(1,2);
% randomly pick half of the other_races with fixed rng
if strcmp(block_names,'E20-1163_0052')
    rng('default');
    
    half_index_black = randsample(size(plot_data{2},1),round(size(plot_data{2},1)/2));
    half_index_white = randsample(size(plot_data{3},1),round(size(plot_data{3},1)/2));
    plot_data_SRvsOR{1} = plot_data{1};
    plot_data_SRvsOR{2} = [plot_data{2}(half_index_black,:);plot_data{3}(half_index_white,:)];
    stats_data_SRvsOR{1} = stats_data{1};
    stats_data_SRvsOR{2} = [stats_data{2}(half_index_black,:);stats_data{3}(half_index_white,:)];
    

%     plot_data_SRvsOR{1} = plot_data{1};
%     plot_data_SRvsOR{2} = [plot_data{2};plot_data{3}];
%     stats_data_SRvsOR{1} = stats_data{1};
%     stats_data_SRvsOR{2} = [stats_data{2};stats_data{3}];
elseif strcmp(block_names,'E20-1163_0106')
    rng('default');
    
    half_index_asian = randsample(size(plot_data{1},1),round(size(plot_data{1},1)/2));
    half_index_black = randsample(size(plot_data{2},1),round(size(plot_data{2},1)/2));
    plot_data_SRvsOR{1} = plot_data{3};
    plot_data_SRvsOR{2} = [plot_data{1}(half_index_asian,:);plot_data{2}(half_index_black,:)];
    stats_data_SRvsOR{1} = stats_data{3};
    stats_data_SRvsOR{2} = [stats_data{1}(half_index_asian,:);stats_data{2}(half_index_black,:)];
    
%     plot_data_SRvsOR{1} = plot_data{3};
%     plot_data_SRvsOR{2} = [plot_data{1};plot_data{2}];
%     stats_data_SRvsOR{1} = stats_data{3};
%     stats_data_SRvsOR{2} = [stats_data{1};stats_data{2}];
else
end
conditions = {'own_race','other_race'}; column = 'condNames9';
%% plot figure based on aboving data
clear h
load('cdcol.mat')
% figureDim = [100 100 .23 .35 ];
figureDim = [1 1 0.36 0.48];%default [0 0 .3 .4]
figure('units', 'normalized', 'outerposition', figureDim)
for ci = 1:length(conditions)
    lineprops.col{1} = plot_params.col(ci,:);
    if ~strcmp(plot_params.eb,'none')
        lineprops.style= '-';
        lineprops.width = plot_params.lw;
        lineprops.edgestyle = '-';
        if strcmp(plot_params.eb,'ste')
            mseb(data_all.time,nanmean(plot_data_SRvsOR{ci}),nanstd(plot_data_SRvsOR{ci})/sqrt(size(plot_data_SRvsOR{ci},1)),lineprops,1);
        
            hold on
        else %'std'
            mseb(data_all.time,nanmean(plot_data_SRvsOR{ci}),nanstd(plot_data_SRvsOR{ci}),lineprops,1);
            hold on
        end
    end
 
    h(ci)=plot(data_all.time,nanmean(plot_data_SRvsOR{ci}),'LineWidth',plot_params.lw,'Color',plot_params.col(ci,:));
    hold on
end
xlim(plot_params.xlim)


if contains(sbj_names{1},'C17_20') && contains(sbj_names{end},'S20_152_HT')
    xlabel(plot_params.xlabel);
else
end


% if strcmp(anat_name,'INSULA')
%     ylabel(plot_params.ylabel);
% else
% end

set(gca,'fontsize',plot_params.textsize)
box off
    


    y_lim = ylim;

    % this part is clusterd based permuation
   
if plot_params.clust_per
    
    ylim_stac = y_lim;  
    indx_per = find(abs(data_all.time-plot_params.clust_per_win(1))<0.001):find(abs(data_all.time-plot_params.clust_per_win(2))<0.001);
    data.time_stac = data_all.time(indx_per);
    
    data_SR = plot_data_SRvsOR{1}(:,indx_per);
    data_SR = permute(data_SR,[3 2 1]);
    data_OR = plot_data_SRvsOR{2}(:,indx_per);
    data_OR  = permute(data_OR,[3 2 1]);

    rng('default')
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_SR, data_OR,[]);   
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of SR higher than OR are ' num2str(clust_info.pos_clust_pval)])
    h(3) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',4,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',20);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of SR lower than OR are ' num2str(clust_info.neg_clust_pval)])
    h(3) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',4,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',20);
    
    
    
   
    
    
    cond_names_stac = cond_names;
    cond_names_stac{1,1} = 'SR';
    cond_names_stac{1,2} = 'OR';
    cond_names_stac{1,3} = 'SR vs OR';
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
leg = legend(h,cond_names_stac,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')


set(gca, 'box', 'off')

%legend off
% set(gca,'XLabel','Time(S)');%chao
xlabel('Time(S)') 

% sites_num = sum(cellfun(@numel, T3{:,'anat'} ));
% sbj_names_num = size(T3,1);
% title([num2str(sites_num),' sites in ' anat_name ' from ',num2str(sbj_names_num),' Subjects'])
