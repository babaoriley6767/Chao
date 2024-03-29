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
side = 'none';%'L','R','none'
site_pick = 'all sites per subj';%all sites per subj

anat_all = {{'INSULA'};{'IFS','IFG'};{'SFS','SFG'}};
anat_name_all = {'INSULA';'IFS IFG';'SFS SFG'};


data_asian_ROL_onset = cell(length(anat_all),1); 
data_black_ROL_onset = cell(length(anat_all),1); 
data_white_ROL_onset = cell(length(anat_all),1); 

cdcol  =[0.9647,0.1647,0; 
         0.3020,0.3922,0.5529;  
         0,0.5800,0.7920]; 


for ai = 1:length(anat_all)
   anat = anat_all{ai};
   anat_name = anat_name_all{ai};
   
   [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');   
   
   switch anat_name
    case 'INSULA'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,2,5,6],:);
            T3.anat{3} = T3.anat{3}(3);
        else
        end
    case 'IFS IFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([5,10],:);
            T3.anat{1} = T3.anat{1}(5);
        else
        end
    case 'SFS SFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,4,6],:);
            T3.anat{2} = T3.anat{2}(3);
            T3.anat{3} = T3.anat{3}(1);
        else
        end
end
disp(channame);

load('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/mat_files/group_diff_encoding_HFB_reall_behav_sheet.mat')
channame_sheet = cell(length(HFB_behav),1);
for i = 1:length(HFB_behav)
    channame_sheet{i} = HFB_behav{i}.channame{1};
end
HFB_behav = HFB_behav(1,ismember(channame_sheet,channame));
HFB_behav_stats = vertcat(HFB_behav{:});

%remove outlier of ROL_onset (modified z-score)

zscorethresh = 3;

idx_clean = HFB_behav_stats.clean_trial==1;


idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);


idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);


ROL_asian = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_asian);
ROL_black = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_black);
ROL_white = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_white);


mZ_column = nan(height(HFB_behav_stats),1);%% zscore
for i = 1:height(HFB_behav_stats)
    if HFB_behav_stats.test_race(i) == 1 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_asian)) ./ mad(ROL_asian,1);
    elseif  HFB_behav_stats.test_race(i) == 2 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_black)) ./ mad(ROL_black,1);
    elseif  HFB_behav_stats.test_race(i) == 3 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_white)) ./ mad(ROL_white,1);
    else
    end
end

mZ_indx = ~(mZ_column>zscorethresh);% key parameter


data_asian_ROL_onset{ai} = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian&mZ_indx,:);

data_black_ROL_onset{ai} = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black&mZ_indx,:);

data_white_ROL_onset{ai} = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white&mZ_indx,:);


end


%%distribution of ROL in different region* condition
figure('Position', [0 0 1100 300])

subplot(1,3,1)
for ai = 1:length(anat_all)
h(ai) = histogram(data_asian_ROL_onset{ai});
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.02;
hold on
end
title('Asian')
set(gca,'fontsize',28)

subplot(1,3,2)
for ai = 1:length(anat_all)
h(ai) = histogram(data_black_ROL_onset{ai});
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.02;
hold on
end
title('Black')
set(gca,'fontsize',28)

subplot(1,3,3)
for ai = 1:length(anat_all)
h(ai) = histogram(data_white_ROL_onset{ai});
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.02;
hold on
end
title('white')
set(gca,'fontsize',28)




data_insula= [data_asian_ROL_onset{1}];
data_vlPFC = [data_asian_ROL_onset{2}];
data_vdPFC = [data_asian_ROL_onset{3};];
data_anova = [data_insula;data_vlPFC;data_vdPFC]';
group_insula = repmat({'insula'},size(data_insula,1),1);
group_vlPFC = repmat({'vlPFC'},size(data_vlPFC,1),1);
group_vdPFC = repmat({'vdPFC'},size(data_vdPFC,1),1);
group = [group_insula;group_vlPFC;group_vdPFC];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_insula),nanstd(data_vlPFC),nanstd(data_vdPFC)];
m1 = multcompare(stats1,'ctype','lsd')

data_insula= [data_black_ROL_onset{1}];
data_vlPFC = [data_black_ROL_onset{2}];
data_vdPFC = [data_black_ROL_onset{3}];
data_anova = [data_insula;data_vlPFC;data_vdPFC]';
group_insula = repmat({'insula'},size(data_insula,1),1);
group_vlPFC = repmat({'vlPFC'},size(data_vlPFC,1),1);
group_vdPFC = repmat({'vdPFC'},size(data_vdPFC,1),1);
group = [group_insula;group_vlPFC;group_vdPFC];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_insula),nanstd(data_vlPFC),nanstd(data_vdPFC)];
m1 = multcompare(stats1,'ctype','lsd')

data_insula= [data_white_ROL_onset{1}];
data_vlPFC = [data_white_ROL_onset{2}];
data_vdPFC = [data_white_ROL_onset{3}];
data_anova = [data_insula;data_vlPFC;data_vdPFC]';
group_insula = repmat({'insula'},size(data_insula,1),1);
group_vlPFC = repmat({'vlPFC'},size(data_vlPFC,1),1);
group_vdPFC = repmat({'vdPFC'},size(data_vdPFC,1),1);
group = [group_insula;group_vlPFC;group_vdPFC];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_insula),nanstd(data_vlPFC),nanstd(data_vdPFC)];
m1 = multcompare(stats1,'ctype','lsd')



%%distribution of ROL in different region* all conditions
values = cell(1,length(anat_all));
edges = cell(1,length(anat_all));
centers = cell(1,length(anat_all));

figure('Position', [0 0 700 300])
subplot(1,2,1)
for ai = 1:length(anat_all)
h(ai) = histogram([data_asian_ROL_onset{ai};data_black_ROL_onset{ai};data_white_ROL_onset{ai}]);
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.03;
h(ai).FaceColor = cdcol(ai,:);
hold on
end
ylabel('Probability')
xlabel('ROL(s)')
leg = legend(h,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',12, 'Interpreter', 'none')
title('all conditions')

subplot(1,2,2)
for ai = 1:length(anat_all)
[values{ai}, edges{ai}] = histcounts([data_asian_ROL_onset{ai};data_black_ROL_onset{ai};data_white_ROL_onset{ai}], 'Normalization', 'probability','Binwidth',.03);
centers{ai} = (edges{ai}(1:end-1)+edges{ai}(2:end))/2;
c(ai) = plot(centers{ai}, values{ai}, '-','color',cdcol(ai,:),'linewidth',1)
hold on
end
xlabel('ROL(s)')
leg = legend(c,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')
title('all conditions')


data_insula= [data_asian_ROL_onset{1};data_black_ROL_onset{1};data_white_ROL_onset{1}];
data_vlPFC = [data_asian_ROL_onset{2};data_black_ROL_onset{2};data_white_ROL_onset{2}];
data_vdPFC = [data_asian_ROL_onset{3};data_black_ROL_onset{3};data_white_ROL_onset{3}];
data_anova = [data_insula;data_vlPFC;data_vdPFC]';
group_insula = repmat({'insula'},size(data_insula,1),1);
group_vlPFC = repmat({'vlPFC'},size(data_vlPFC,1),1);
group_vdPFC = repmat({'vdPFC'},size(data_vdPFC,1),1);
group = [group_insula;group_vlPFC;group_vdPFC];
%cd('/Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/globe_analysis_figures');%plz adjust accordingly
[p1,tbl1,stats1] = anova1(data_anova,group);
std1 = [nanstd(data_insula),nanstd(data_vlPFC),nanstd(data_vdPFC)];
m1 = multcompare(stats1,'ctype','lsd')




figure('Position', [0 0 700 300])
subplot(1,2,1)
for ai = 1:length(anat_all)
h(ai) = histogram([data_asian_ROL_onset{ai}]);
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.03;
h(ai).FaceColor = cdcol(ai,:);
hold on
end
ylabel('Probability')
xlabel('ROL(s)')
leg = legend(h,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',12, 'Interpreter', 'none')
title('Asian')

subplot(1,2,2)
for ai = 1:length(anat_all)
[values{ai}, edges{ai}] = histcounts([data_asian_ROL_onset{ai}], 'Normalization', 'probability','Binwidth',.03);
centers{ai} = (edges{ai}(1:end-1)+edges{ai}(2:end))/2;
c(ai) = plot(centers{ai}, values{ai}, '-','color',cdcol(ai,:),'linewidth',1)
hold on
end
xlabel('ROL(s)')
leg = legend(c,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')
title('Asian')



figure('Position', [0 0 700 300])
subplot(1,2,1)
for ai = 1:length(anat_all)
h(ai) = histogram([data_black_ROL_onset{ai}]);
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.03;
h(ai).FaceColor = cdcol(ai,:);
hold on
end
ylabel('Probability')
xlabel('ROL(s)')
leg = legend(h,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',12, 'Interpreter', 'none')
title('Black')

subplot(1,2,2)
for ai = 1:length(anat_all)
[values{ai}, edges{ai}] = histcounts([data_black_ROL_onset{ai}], 'Normalization', 'probability','Binwidth',.03);
centers{ai} = (edges{ai}(1:end-1)+edges{ai}(2:end))/2;
c(ai) = plot(centers{ai}, values{ai}, '-','color',cdcol(ai,:),'linewidth',1)
hold on
end
xlabel('ROL(s)')
leg = legend(c,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')
title('Black')



figure('Position', [0 0 700 300])
subplot(1,2,1)
for ai = 1:length(anat_all)
h(ai) = histogram([data_white_ROL_onset{ai}]);
h(ai).Normalization = 'probability';
h(ai).BinWidth = 0.03;
h(ai).FaceColor = cdcol(ai,:);
hold on
end
ylabel('Probability')
xlabel('ROL(s)')
leg = legend(h,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',12, 'Interpreter', 'none')
title('White')

subplot(1,2,2)
for ai = 1:length(anat_all)
[values{ai}, edges{ai}] = histcounts([data_white_ROL_onset{ai}], 'Normalization', 'probability','Binwidth',.03);
centers{ai} = (edges{ai}(1:end-1)+edges{ai}(2:end))/2;
c(ai) = plot(centers{ai}, values{ai}, '-','color',cdcol(ai,:),'linewidth',1)
hold on
end
xlabel('ROL(s)')
leg = legend(c,anat_name_all,'Location','Northeast', 'AutoUpdate','off');%cond_names has the trial infomation(default), and cond_names2 is about the category
legend boxoff
set(leg,'fontsize',plot_params.legendfontsize, 'Interpreter', 'none')
title('White')
% title('Asian')
% set(gca,'fontsize',28)