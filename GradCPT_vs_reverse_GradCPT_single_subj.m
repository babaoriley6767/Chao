%get ready with the addpath,plz adjust accordingly
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

project_name ='GradCPT';% 'race_encoding_simple'or'Grad_CPT'  % The might be error here because of the naming some are GradCPT and some are Grad_CPT
project_name_r = 'GradCPT_re';
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = false;% clusterd based permuation
plot_params.clust_per_win = [0 1];
%%
%make a specific selection of conditions and coloring

conditions = {'mtn','city'};column = 'condNames';
load cdcol.mat
plot_params.col = [cdcol.carmine;cdcol.ultramarine];
%%
%concatenate  data for conditions
plot_data = cell(1,length(conditions));
plot_data_r = cell(1,length(conditions));
stats_data = cell(1,length(conditions));
stats_data_r = cell(1,length(conditions));

sbj_name = 'S20_158';%set basic pipeline parameters
if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
    center = 'China';
else
    center = 'Stanford';
end
block_names = BlockBySubj(sbj_name,project_name);
block_names = block_names(1,5:8)
block_names_r = BlockBySubj(sbj_name,project_name_r)
dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
dirs_r = InitializeDirs(project_name_r,sbj_name,comp_root,server_root,code_root);
load([dirs.data_root,filesep,'OriginalData',filesep,sbj_name,filesep,'global_',project_name,'_',sbj_name,'_',block_names{1},'.mat'])
for i = 1:length(globalVar.channame)
    
    
    data_all = concatBlocks(sbj_name, block_names,dirs,[81],'HFB','Band',{'wave'},['stimlock']);%'stimlock_bl_corr'
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
        plot_data{ci} = data_all.wave_sm(grouped_trials{ci},:);% we use smoothed data for plotting
        stats_data{ci} = data_all.wave(grouped_trials{ci},:); % this part of data were prepared for further comparison (original non-smoothed data)
    end
    
    
    
    data_all_r = concatBlocks(sbj_name, block_names_r,dirs_r,[81],'HFB','Band',{'wave'},['stimlock']);%'stimlock_bl_corr'
    %smooth is in the trial level
    winSize = floor(data_all_r.fsample*plot_params.sm);
    gusWin = gausswin(winSize)/sum(gausswin(winSize));
    data_all_r.wave_sm = convn(data_all_r.wave,gusWin','same');
    [grouped_trials_all,~] = groupConds(conditions,data_all_r.trialinfo,column,'none',[],false);
    [grouped_trials,cond_names] = groupConds(conditions,data_all_r.trialinfo,column,stats_params.noise_method,stats_params.noise_fields_trials,false);
    % this part is to exclude HFB over a fixed threshold
    if plot_params.single_trial_replot
        thr_raw =[];
        for di = 1:size(data_all_r.wave,1)
            if ~isempty(find(data_all_r.wave(di,:)>=plot_params.single_trial_thr));
                fprintf('You have deleted the data over threshold %d from the data \n',plot_params.single_trial_thr);
            else
            end
        end
        [thr_raw,thr_column] = find(data_all_r.wave >= plot_params.single_trial_thr);
        thr_raw = unique(thr_raw);
    end
    for ci = 1:length(conditions)
        grouped_trials{ci} = setdiff(grouped_trials{ci},thr_raw);% make the grouped_trial and thr_raw in together
        plot_data_r{ci} = data_all_r.wave_sm(grouped_trials{ci},:);% we use smoothed data for plotting
        stats_data_r{ci} = data_all_r.wave(grouped_trials{ci},:); % this part of data were prepared for further comparison (original non-smoothed data)
    end
    
    
    
    %build data   
    plot_data_syn = cell(1,2);
    stats_data_syn = cell(1,2);
    plot_data_syn{1,1} = plot_data{1,1};
    plot_data_syn{1,2} = plot_data_r{1,1};
    stats_data_syn{1,1} = stats_data{1,1};
    stats_data_syn{1,2} = stats_data_r{1,1};
    
  
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
            mseb(data_all.time,nanmean(plot_data_syn{ci}),nanstd(plot_data_syn{ci})/sqrt(size(plot_data_syn{ci},1)),lineprops,1);
            
            hold on
        else %'std'
            mseb(data_all.time,nanmean(plot_data_syn{ci}),nanstd(plot_data_syn{ci}),lineprops,1);
            hold on
        end
    end
    
    h(ci)=plot(data_all.time,nanmean(plot_data_syn{ci}),'LineWidth',plot_params.lw,'Color',plot_params.col(ci,:));
    hold on
end
xlim(plot_params.xlim)



set(gca,'fontsize',plot_params.textsize)
box off


y_lim = ylim;

% this part is clusterd based permuation

if plot_params.clust_per
    
    ylim_stac = y_lim;
    indx_per = find(abs(data_all.time-plot_params.clust_per_win(1))<0.001):find(abs(data_all.time-plot_params.clust_per_win(2))<0.001);
    data.time_stac = data_all.time(indx_per);
    
    data_mtn = stats_data{1}(:,indx_per);
    data_mtn = permute(data_mtn,[3 2 1]);
    data_city = stats_data{2}(:,indx_per);
    data_city = permute(data_city,[3 2 1]);
    rng('default')
    
    [pval, t_orig, clust_info, seed_state, est_alpha] = clust_perm2(data_mtn, data_city,[]);
    
    cluster_indx = find(clust_info.pos_clust_pval<0.05);
    cluster_sig = ismember(clust_info.pos_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of mtn higher than city are ' num2str(clust_info.pos_clust_pval)])
    h(3) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    cluster_indx = find(clust_info.neg_clust_pval<0.05);
    cluster_sig = ismember(clust_info.neg_clust_ids,cluster_indx)*1;
    cluster_sig(cluster_sig==0)=NaN;
    warning(['cluster p value of mtn lower than city are ' num2str(clust_info.neg_clust_pval)])
    h(3) = plot(data.time_stac, cluster_sig.*ylim_stac(2)-0.05, '-*','Color',cdcol.black,'LineWidth',2,'MarkerIndices',1:35:length(data.time_stac),'MarkerSize',14);
    
    cond_names_stac = cond_names;
    cond_names_stac{1,1} = 'mtn';
    cond_names_stac{1,2} = 'city';
    cond_names_stac{1,3} = 'mtn vs city';
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
leg = legend(h,{'GradCPT_mtn','reverse_GradCPT_mtn'},'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')

%legend off
% set(gca,'XLabel','Time(S)');%chao
xlabel('Time(S)')

%get rid of the upper and right axis
set(gca, 'box', 'off')



suptitle([data_all.label])

dir_out = '/Users/chao/Desktop/158-GradCPTvsRE';
folder_name = 'GradCPTvsRE';
fn_out = sprintf('%s/%s_%s_%s.png',dir_out,sbj_name,data_all.label,folder_name);%
savePNG(gcf, 200, fn_out)

close
end


