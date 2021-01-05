
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
            idx2 = T{i}.group_diff;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation/T{i}.group_diff %default
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
            idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
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
loc=cellfun('isempty', T3{:,'anat'} );
T3(loc,:)=[];
% % s = rng;
% rand37 = randperm(24,7);
% T3.anat{7} = T3.anat{7}(rand37);
%%
%define the plot and stats parameters first
project_name ='race_encoding_simple';% 'race_encoding_simple'or'Grad_CPT'
plot_params = genPlotParams(project_name,'timecourse');
plot_params.single_trial_replot = true;
plot_params.single_trial_thr = 15;%the threshold of HFB it could be like 10 15 20 ...
stats_params = genStatsParams(project_name);
plot_params.single_trial = false;
plot_params.clust_per = true;% clusterd based permuation
%%
%make a specific selection of conditions and coloring
if strcmp(project_name,'Grad_CPT')
    conditions = {'mtn','city'};column = 'condNames';
    load cdcol.mat
    plot_params.col = [cdcol.carmine;cdcol.ultramarine];
elseif strcmp(project_name,'race_encoding_simple')
    conditions = {'asian','black','white'}; column = 'condNames';
   % conditions = {'own_race','other_races'};column = 'condNames9';
    %conditions = {'my_race_ans','not_my_race_ans'};column = 'condNames8';
%     load cdcol.mat
%     plot_params.col = [cdcol.carmine;cdcol.ultramarine];
end

%%
%concatenate  data for conditions

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
            
            if i == 1
                count_indx = j;
            elseif i > 1
                count_indx = sum(cellfun(@numel,T3.anat(1:(i-1),:)))+j;
            else
            end
            
            
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
                    stats_data{ci} = [stats_data{ci};data_all.wave(grouped_trials{ci},:)]; % this part of data were prepared for further comparison (original non-smoothed data)
                    stats_data_all{ci} = [stats_data_all{ci};data_all.wave(grouped_trials_all{ci},:)];
                    cond = conditions{ci};
                    stats_data_ROL.(cond){count_indx,1}=data_all.wave(grouped_trials{ci},:);%this is the ROL data
                end
                
        end
    else
    end
end


% sum(cellfun(@numel,T3.anat))
% % randomly pick half of the other_races with fixed rng
% if strcmp(column,'condNames9')
%     rng('default');
%     half_index = randsample(size(plot_data{2},1),round(size(plot_data{2},1)/2));
%     plot_data{2} = plot_data{2}(half_index,:);
% else
% en

%% ROL from Jessica and Omri NC parameters
ROL_params = genROLParams_NC(project_name)% set the parameters
%getROL break down
nstim = 1;
fs = data_all.fsample; % sampling rate
stimtime = 0;%%%%%% time 0
stiminds = find(data_all.time>stimtime,1); %%%%%%% indx of time 0
count_all = sum(cellfun(@numel,T3.anat));%%%% number of total sites in the ROL analysis
befInd = round(ROL_params.pre_event * fs);%%% indx numbers before time 0
aftInd = round(ROL_params.dur * fs);%%%% indx numbers after time 0
time = (befInd:aftInd)/fs;%% time series
%% calculate ROL across sites
for ei = 1:count_all
    for ci = 1:length(conditions)
        cond = conditions{ci};
         stats_data_ROL.(cond){ei} = stats_data_ROL.(cond){ei} (:,stiminds+befInd:stiminds+aftInd);%% define the ROL data length and range
    end
end


for ci = 1:length(conditions)% data initiation
    cond = conditions{ci};
    ROL.(cond).peaks = cell(count_all,nstim);
    ROL.(cond).onsets = cell(count_all,nstim);
end

% ROL core calc
for ei = 1:sum(cellfun(@numel,T3.anat))% ROL main part
    for ci = 1:length(conditions)
        cond = conditions{ci};
        for ii = 1:nstim
            data.wave = stats_data_ROL.(cond){ei,ii};
            data.time = time;
            [Resp_data]= ROLbootstrap_NC(data, ROL_params);
            ROL.(cond).peaks{ei,ii} = Resp_data.peaks;
            ROL.(cond).onsets{ei,ii} = Resp_data.onsets;
        end
    end
    disp(['Computing ROL for sites: ',num2str(ei)])
end
%% calculate ROL across all trials

    for ci = 1:length(conditions)
         stats_data{ci}  = stats_data{ci}(:,stiminds+befInd:stiminds+aftInd);%% define the ROL data length and range
    end
    
        for ci = 1:length(conditions)
            cond = conditions{ci};
            data.wave = stats_data{ci};
            data.time = time;
            [Resp_data]= ROLbootstrap_NC(data, ROL_params);
            ROL_all_trials.(cond).peaks = Resp_data.peaks;
            ROL_all_trials.(cond).onsets = Resp_data.onsets;
            disp(['Computing ROL_all_trials for condition: ',num2str(ci)])
        end


%%
%Sites for which a ROL value could not be obtained in 50% of the trials or more were discarded from the analysis.
% range first 50% second place

%get rid of the sites with too much nan and make a index of nan?the
%following is based on each condition
idx_nan = ones(length(channame),length(conditions));
for elec_inspect = 1: length(channame)
    for ci = 1:length(conditions)
        if sum(isnan(ROL.(conditions{ci}).onsets{elec_inspect}))>round(length(ROL.(conditions{ci}).onsets{elec_inspect}))*.5
            idx_nan(elec_inspect,ci) = 0;
        else
        end
    end
end
idx_nan = cumprod(idx_nan,2);
idx_nan = idx_nan(:,end);
idx_nan = logical(idx_nan)

%get rid of the sites with too much nan and make a index of nan?the
%following is based on each electrode
idx_nan = ones(length(channame),1);
for elec_inspect = 1: length(channame)
        if sum(isnan([ROL.asian.onsets{elec_inspect} ROL.black.onsets{elec_inspect} ROL.white.onsets{elec_inspect}]))>round(length([ROL.asian.onsets{elec_inspect} ROL.black.onsets{elec_inspect} ROL.white.onsets{elec_inspect}]))*.5
            idx_nan(elec_inspect) = 0;
        else
        end
end
idx_nan = logical(idx_nan)

%take care of the time rage that of interest
twin = [0.1 1];%[0.1 1]
for elec_inspect = 1: length(channame)
    for ci = 1:length(conditions)
        ROL.(conditions{ci}).onsets{elec_inspect}(ROL.(conditions{ci}).onsets{elec_inspect} < twin(1) | ROL.(conditions{ci}).onsets{elec_inspect} > twin(2)) = nan;
    end
end


%creat little version of ROL and channel names
for ci = 1:length(conditions)
ROL.(conditions{ci}).onsets_lit = ROL.(conditions{ci}).onsets(idx_nan,:);
ROL.(conditions{ci}).peaks_lit = ROL.(conditions{ci}).onsets(idx_nan,:);
end
channame_lit = channame(idx_nan,:);
ROL.channame = channame;
ROL.channame_lit = channame_lit;

%%
% plot each site
for elec_inspect = 1: length(channame_lit)
    %chao plot the distribution of ROL.thr
%     figure(elec_inspect)
    figure('units', 'normalized', 'outerposition', [0 0 0.2 0.4])
    for ci = 1:length(conditions)
        h(ci) = plot(sort(ROL.(conditions{ci}).onsets_lit{elec_inspect},'descend'), 1:length(ROL.(conditions{ci}).onsets_lit{elec_inspect}),'o','Color',plot_params.col(ci,:),'LineWidth',1);
        hold on
    end
    title('Single trial ROL NC')
    suptitle(channame_lit{elec_inspect})
    xlabel('ROL (s)')
    xlim([0 1])
    xticks(0:0.1:2)
    xticklabels(0:0.1:2)    
    ylabel('Trials')
    set(gca,'fontsize',10)
    grid on
    leg = legend(h,conditions,'Location','Northeast', 'AutoUpdate','off');%
    legend boxoff%
    set(leg,'fontsize',18, 'Interpreter', 'none')
end
%%
%anova and post hoc
stats.p = cell(length(channame_lit),1);
stats.tbl = cell(length(channame_lit),1);
stats.stats = cell(length(channame_lit),1);
stats.m = cell(length(channame_lit),1);
stats.std = cell(length(channame_lit),1);

% stats for each sites
for elec_inspect = 1: length(channame_lit)
data_asian = ROL.asian.onsets_lit{elec_inspect}';
data_black = ROL.black.onsets_lit{elec_inspect}';
data_white = ROL.white.onsets_lit{elec_inspect}';
data_anova = [data_asian;data_black;data_white];
group_asian = repmat({'asian'},size(data_asian,1),1);
group_black = repmat({'black'},size(data_black,1),1);
group_white = repmat({'white'},size(data_white,1),1);
group = [group_asian;group_black;group_white];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[stats.p{elec_inspect},stats.tbl{elec_inspect},stats.stats{elec_inspect}] = anova1(data_anova,group);
stats.std{elec_inspect} = [nanstd(data_asian),nanstd(data_black ),nanstd(data_white)];
stats.m{elec_inspect} = multcompare(stats.stats{elec_inspect},'ctype','lsd');
% m1 = multcompare(stats1,'ctype','bonferroni')
end

% stats across sites
data_asian = cellfun(@nanmean,ROL.asian.onsets_lit);
data_black = cellfun(@nanmean,ROL.black.onsets_lit);
data_white = cellfun(@nanmean,ROL.white.onsets_lit);
data_anova = [data_asian;data_black;data_white];
group_asian = repmat({'asian'},size(data_asian,1),1);
group_black = repmat({'black'},size(data_black,1),1);
group_white = repmat({'white'},size(data_white,1),1);
group = [group_asian;group_black;group_white];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_asian),nanstd(data_black ),nanstd(data_white)];
m1 = multcompare(stats1,'ctype','lsd')
% m1 = multcompare(stats1,'ctype','bonferroni')

%stats across trials
data_asian = horzcat(ROL.asian.onsets_lit{:,:})';
data_asian = data_asian(~isnan(data_asian));

data_black = horzcat(ROL.black.onsets_lit{:,:})';
data_black = data_black(~isnan(data_black));

data_white = horzcat(ROL.white.onsets_lit{:,:})';
data_white = data_white(~isnan(data_white));

data_anova = [data_asian;data_black;data_white]';
group_asian = repmat({'asian'},size(data_asian,1),1);
group_black = repmat({'black'},size(data_black,1),1);
group_white = repmat({'white'},size(data_white,1),1);
group = [group_asian;group_black;group_white];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_asian),nanstd(data_black ),nanstd(data_white)];
m1 = multcompare(stats1,'ctype','lsd')


%%
save('Race_ROL_group_diff_dlPFC.mat','ROL');
%%
fn_out = sprintf('%s/%s_%s_%s_ROL2_fig.png',dir_out,sbj_name,subjVar.labels_EDF{elec_inspect},project_name);%
savePNG(gcf, 100, fn_out)%%%%%%%%%

% Exlucde (set to nan) ROLS beyond the temporal window of choice and calculate mean and median ROL across trials per electrode
for i = 1:length(ROL_NC_fig.(sub_cond{1}).onsets)%length of the elecs
    for j = 1:length(sub_cond);
    ROL_NC_fig.(sub_cond{j}).onsets{i}(ROL_NC_fig.(sub_cond{j}).onsets{i} < twin(1) | ROL_NC_fig.(sub_cond{j}).onsets{i} > twin(2)) = nan;
    ROL_NC_fig.(sub_cond{j}).median(i) = nanmedian(ROL_NC_fig.(sub_cond{j}).onsets{i});
    ROL_NC_fig.(sub_cond{j}).sd(i) = nanstd(ROL_NC_fig.(sub_cond{j}).onsets{i});
    ROL_NC_fig.(sub_cond{j}).mean(i) = nanmean(ROL_NC_fig.(sub_cond{j}).onsets{i});
    end
end

% save the new ROL as ROL_fig in the result folder

fn_out = sprintf('%s%s_%s_ROL_NC_fig2.mat',dir_out,sbj_name,project_name);%%%%%%%%%
save(fn_out,'ROL_NC_fig');

for i = 1:length(sub_cond)
    sprintf('%s_%s_%s_ROL_NC_mean_of_%s_is_%d',sbj_name,subjVar.labels_EDF{elec_inspect},project_name,sub_cond{i},(ROL_NC_fig.(sub_cond{i}).mean(elec_inspect)*1000))
end

for i = 1:length(sub_cond)
    sprintf('%s_%s_%s_ROL_NC_ratio_of_%s_from_%d_trials_is_%d',sbj_name,subjVar.labels_EDF{elec_inspect},project_name,sub_cond{i},length(ROL_NC_fig.(sub_cond{i}).onsets{elec_inspect,1}),sum(~isnan(ROL_NC_fig.(sub_cond{i}).onsets{elec_inspect,1}))*100/length(ROL_NC_fig.(sub_cond{i}).onsets{elec_inspect,1}) )
end

patlist = {'S1_X7R','S1_X6R','S1_X5R','S2_X6R','S3_X6R','S4_X7L'};
numlist = [72 71 70 6 176 86];
i = find(abs(numlist-elec_inspect)<0.001);

asian{i} = ROL_NC_fig.asian.onsets{numlist(i)};
black{i} = ROL_NC_fig.black.onsets{numlist(i)};
white{i} = ROL_NC_fig.white.onsets{numlist(i)};
ind1 = ~isnan(asian{i});
ind2 = ~isnan(black{i});
ind3 = ~isnan(white{i});
asian{i} = asian{i}(ind1)';
black{i} = black{i}(ind2)';
white{i} = white{i}(ind3)';
lna = length(asian{i});
lnb = length(black{i});
lnw = length(white{i});
hfb{i} = [asian{i};black{i};white{i}];
sites{i} = cell(lna+lnb+lnw,1);
sites{i}(:) = {patlist{i}};%% tricky
race1 = cell(lna,1);
race1(:) = {'asian'};
race2 = cell(lnb,1);
race2(:) = {'black'};
race3 = cell(lnw,1);
race3(:) = {'white'};
race{i} = vertcat(race1,race2,race3);

hfbcat = vertcat(hfb{:});
racecat = vertcat(race{:});
sitescat = vertcat(sites{:});

asian_mean = nanmean(data{1});
asian_sd = nanstd(data{1});
black_mean = nanmean(data{2});
black_sd = nanstd(data{2});
white_mean = nanmean(data{3});
white_sd = nanstd(data{3});


for i = 1:6
[h1,p1] = ttest2(hfb{i}(strcmp(race{i},'asian')),hfb{i}(strcmp(race{i},'black')));
[h2,p2] = ttest2(hfb{i}(strcmp(race{i},'asian')),hfb{i}(strcmp(race{i},'white')));
[h3,p3] = ttest2(hfb{i}(strcmp(race{i},'black')),hfb{i}(strcmp(race{i},'white')));
end

[h1,p1] = ttest2(hfbcat(strcmp(racecat,'asian')),hfbcat(strcmp(racecat,'black')));
[h2,p2] = ttest2(hfbcat(strcmp(racecat,'asian')),hfbcat(strcmp(racecat,'white')));
[h3,p3] = ttest2(hfbcat(strcmp(racecat,'black')),hfbcat(strcmp(racecat,'white')));

save('ROL_group_race.mat','hfb','sites','race','hfbcat','racecat','sitescat');
table_ROL = table(racecat,sitescat,hfbcat)
writetable(table_ROL,'table_ROL.xlsx','Sheet','1')

asian_any = vertcat(asian{1},asian{2},asian{3},asian{4},asian{5},asian{6});
black_any = vertcat(black{1},black{2},black{3},black{4},black{5},black{6});
white_any = vertcat(white{1},white{2},white{3},white{4},white{5},white{6});
asian_race = vertcat(asian{1},asian{4},asian{5},asian{6});
black_race = vertcat(black{1},black{4},black{5},black{6});
white_race = vertcat(white{1},white{4},white{5},white{6});
[h1,p1] = ttest2(asian_any,black_any);
[h2,p2] = ttest2(asian_any,white_any);
[h3,p3] = ttest2(black_any,white_any);
[h4,p4] = ttest2(asian_race,black_race);
[h5,p5] = ttest2(asian_race,white_race);
[h6,p6] = ttest2(black_race,white_race);

figure
cats = categorical({'asian','black','white'});
y = [asian_mean;black_mean;white_mean];
sd = [asian_sd,black_sd,white_sd];
hold on
for i = 1:3
b(i) = barh(cats(i),y(i));
b(i).FaceColor = plot_params.col(i,:) ;
er(i) = errorbar(y(i),cats(i),sd(i),sd(i),'.','horizontal');
er(i).Color = [0 0 0];er(i).LineWidth = 3;
end
set(gca,'Xlim',[0 1.5])
set(gca,'FontSize',26)
xlabel('ROL')
% ylabel('RACE')
x0      = 1800;   % Screen position
y0      = 800;   % Screen position
width  = 600; % Width of figure
height = 200; % Height of figure (by default in pixels)
set(gcf,'Position', [x0 y0 width height]);
