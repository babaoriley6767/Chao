%%
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

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
HFB_params.ROL_onsets = [0 1];
HFB_params.ROL_peaks = [.15 .15];




sbj_name = 'S20_158';%set basic pipeline parameters
if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
    center = 'China';
else
    center = 'Stanford';
end
block_names = BlockBySubj(sbj_name,project_name);
block_names = block_names(1,3);%%%%%!!!!!!!!!!!!!!!!!key parameters!!!!!!!!!!!!!!!!!!!
dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);



data_all = concatBlocks(sbj_name, block_names,dirs,[81],'HFB','Band',{'wave'},['stimlock_bl_corr']);%'stimlock_bl_corr'


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


% Exclude values earlier than 10ms
ROL_onsets(ROL_onsets < ROL_params.range(1)) = nan;
ROL_onsets(ROL_onsets > ROL_params.range(2)) = nan;
% ?Sites for which a ROL value could not be obtained in 50% of the trials or more were discarded from the analysis.




data_SRvsOR = cell(1,2);
ROL_onset_SRvsOR = cell(1,2);
time_range = dsearchn(data_all.time',[0 1]');

% randomly pick half of the other_races with fixed rng
if strcmp(block_names,'E20-1163_0052')||strcmp(block_names,'E20-1163_0091')
    rng('default');
    
    half_index_black = sort(randsample(size(grouped_trials{2},1),round(size(grouped_trials{2},1)/2)));
    half_index_black = grouped_trials{2}(half_index_black);
    half_index_white = sort(randsample(size(grouped_trials{3},1),round(size(grouped_trials{3},1)/2)));
    half_index_white = grouped_trials{3}(half_index_white);
    other_race_index = sort([half_index_black;half_index_white]);
    
    other_race_index = sort([grouped_trials{2};grouped_trials{3}]);
    
    data_SRvsOR{1} = data_all.wave(grouped_trials{1},time_range(1):time_range(2));
    data_SRvsOR{2} = data_all.wave(other_race_index,time_range(1):time_range(2));
    
    ROL_onset_SRvsOR{1} = ROL_onsets(grouped_trials{1});
    ROL_onset_SRvsOR{2} = ROL_onsets(other_race_index);
    

%     plot_data_SRvsOR{1} = plot_data{1};
%     plot_data_SRvsOR{2} = [plot_data{2};plot_data{3}];
%     stats_data_SRvsOR{1} = stats_data{1};
%     stats_data_SRvsOR{2} = [stats_data{2};stats_data{3}];
elseif strcmp(block_names,'E20-1163_0106')
    rng('default');
    
    half_index_asian = sort(randsample(size(grouped_trials{1},1),round(size(grouped_trials{1},1)/2)));
    half_index_asian = grouped_trials{1}(half_index_asian);
    half_index_black = sort(randsample(size(grouped_trials{2},1),round(size(grouped_trials{2},1)/2)));
    half_index_black = grouped_trials{2}(half_index_black);
    other_race_index = sort([half_index_asian;half_index_black]);
    
    other_race_index = sort([grouped_trials{1};grouped_trials{2}]);
    
    data_SRvsOR{1} = data_all.wave(grouped_trials{3},time_range(1):time_range(2));
    data_SRvsOR{2} = data_all.wave(other_race_index,time_range(1):time_range(2));
    
    ROL_onset_SRvsOR{1} = ROL_onsets(grouped_trials{3});
    ROL_onset_SRvsOR{2} = ROL_onsets(other_race_index);
    
%     plot_data_SRvsOR{1} = plot_data{3};
%     plot_data_SRvsOR{2} = [plot_data{1};plot_data{2}];
%     stats_data_SRvsOR{1} = stats_data{3};
%     stats_data_SRvsOR{2} = [stats_data{1};stats_data{2}];
else
end


%remove outlier of ROL_onset (modified z-score)

zscorethresh = 3;


ROL_onset_SRvsOR_mZ = ROL_onset_SRvsOR;
mZ = ROL_onset_SRvsOR;
data_SRvsOR_mZ = data_SRvsOR;


for i = 1:length(ROL_onset_SRvsOR)
    for j = 1:length(ROL_onset_SRvsOR{i})
        mZ{i}(j) = norminv(.75)*(ROL_onset_SRvsOR{i}(j)-nanmedian(ROL_onset_SRvsOR{i}))./mad(ROL_onset_SRvsOR{i});
        if mZ{i}(j) > zscorethresh
            ROL_onset_SRvsOR_mZ{i}(j) = nan;
        end
    end
end

%% data converting
for i = 1:length(data_SRvsOR)
    data_SRvsOR_mZ{i} = data_SRvsOR{i}(~isnan(ROL_onset_SRvsOR_mZ{i}),:);
    ROL_onset_SRvsOR_mZ{i} = ROL_onset_SRvsOR{i}(~isnan(ROL_onset_SRvsOR_mZ{i}));
end

%%
% time-trial-hfb plot of ROL onset asian black and white

time_range = dsearchn(data_all.time',[0 1]');


[ROL_SR,SR_sort] = sort(ROL_onset_SRvsOR_mZ{1});
data_SR = data_SRvsOR_mZ{1}(SR_sort,:);
[ROL_OR,OR_sort] = sort(ROL_onset_SRvsOR_mZ{2});
data_OR = data_SRvsOR_mZ{2}(OR_sort,:);


% plot results

figure('Position', [1500 500 450 600]),clf
%  colormap(redblue)
subplot(2,1,1)
contourf(data_all.time(time_range(1):time_range(2)),size(data_SR,1):-1:1,data_SR,40,'linecolor','none')
set(gca,'clim',[-1.5 1.5])
hold on
plot(ROL_SR,size(data_SR,1):-1:1,'k--','LineWidth',4.5)
title('SR, Time-Trial-HFB plot, sorted by ROL')
set(gca,'fontsize',18)
ylabel('trials')
h = colorbar;
ylabel(h, 'HFB Z-score')

subplot(2,1,2)
contourf(data_all.time(time_range(1):time_range(2)),size(data_OR,1):-1:1,data_OR,40,'linecolor','none')
set(gca,'clim',[-1.5 1.5])
hold on
plot(ROL_OR,size(data_OR,1):-1:1,'k--','LineWidth',4.5)
title('OR, Time-Trial-HFB plot, sorted by ROL')
set(gca,'fontsize',18)
ylabel('trials')
xlabel('times(S)')
h = colorbar;
ylabel(h, 'HFB Z-score')


%MannWhitney test
[p,h,stats] = ranksum(ROL_SR,ROL_OR);

%%
data_stats = [ROL_OR;ROL_SR];
group_OR = cell(length(ROL_OR),1);
group_OR(:) = deal({'OR'});
group_SR = cell(length(ROL_SR),1);
group_SR(:) = deal({' SR'});%add space this data column will come to the first
data_groups = vertcat(group_OR,group_SR);
clear g

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups,'linestyle','--');
g(1,1).stat_violin('fill','transparent');
g(1,1).set_title('stat_violin()');

g(1,1)=gramm('x',data_groups,'y',data_stats,'color',data_groups);
g(1,1).set_names('x',[],'y','ROL(S)','color','Origin');
g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.2);
g(1,1).set_color_options('map',plot_params.col);
g(1,1).axe_property('FontSize',20)


figure('Position',[100 100 600 550]);
g.draw();

[h,p,ci,stats] = ttest2(ROL_SR,ROL_OR)

mean(ROL_SR)
std(ROL_SR)

405


