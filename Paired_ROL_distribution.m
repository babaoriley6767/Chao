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

anat_all = {{'INSULA'};{'IFS','IFG'}};
anat_name_all = {'INSULA';'IFS IFG'};

cdcol  =[0.9647,0.1647,0; 
         0.3020,0.3922,0.5529;  
         0,0.5800,0.7920]; 

% data_asian_ROL_onset = cell(length(anat_all),1); 
% data_black_ROL_onset = cell(length(anat_all),1); 
% data_white_ROL_onset = cell(length(anat_all),1); 
% 
% cdcol  =[0.9647,0.1647,0; 
%          0.3020,0.3922,0.5529;  
%          0,0.5800,0.7920]; 

channame = cell(1,length(anat_all));
for ai = 1:length(anat_all)
   anat = anat_all{ai};
   anat_name = anat_name_all{ai};
   [channame{ai},T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');   
   
   switch anat_name
    case 'INSULA'
        if strcmp(site_pick,'One site  per subj')
            channame{ai} = channame([1,2,5,6],:);
            T3.anat{3} = T3.anat{3}(3);
        else
        end
    case 'IFS IFG'
        if strcmp(site_pick,'One site  per subj')
            channame{ai} = channame([5,10],:);
            T3.anat{1} = T3.anat{1}(5);
        else
        end
    case 'SFS SFG'
        if strcmp(site_pick,'One site  per subj')
            channame{ai} = channame([1,4,6],:);
            T3.anat{2} = T3.anat{2}(3);
            T3.anat{3} = T3.anat{3}(1);
        else
        end
end
disp(channame{ai});
end

idx = [];
for ch1 = 1:length(channame{1})
    for ch2 = 1:length(channame{2})
        if contains(channame{1}{ch1},channame{2}{ch2}(1:6))
            idx=[idx;ch1,ch2];
        else
        end
    end
end


%make a big channame sheet
load('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/mat_files/group_diff_encoding_HFB_reall_behav_sheet.mat')
channame_sheet = cell(length(HFB_behav),1);
for i = 1:length(HFB_behav)
    channame_sheet{i} = HFB_behav{i}.channame{1};
end
ROL_all = [];
ROL_asian = [];
ROL_black = [];
ROL_white = [];
for i = 1:length(idx)
channame_pair = {channame{1}{idx(i,1)};channame{2}{idx(i,2)}};
HFB_behav_pair = HFB_behav(1,ismember(channame_sheet,channame_pair));
HFB_behav_stats = vertcat(HFB_behav_pair{:});

idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);
idx_1st = [true(108,1);false(108,1)];
idx_2nd = [false(108,1);true(108,1)];

ROL_all = [ROL_all;HFB_behav_stats.ROL_onsets(idx_1st)-HFB_behav_stats.ROL_onsets(idx_2nd)];
ROL_asian = [ROL_asian;HFB_behav_stats.ROL_onsets(idx_1st&idx_asian)-HFB_behav_stats.ROL_onsets(idx_2nd&idx_asian)];
ROL_black = [ROL_black;HFB_behav_stats.ROL_onsets(idx_1st&idx_black)-HFB_behav_stats.ROL_onsets(idx_2nd&idx_black)];
ROL_white = [ROL_white;HFB_behav_stats.ROL_onsets(idx_1st&idx_white)-HFB_behav_stats.ROL_onsets(idx_2nd&idx_white)];
end
%remove outlier of ROL_onset (modified z-score)
ROL_all = ROL_all(~isnan(ROL_all));

figure('Position', [0 0 400 400])
h = histogram(ROL_all);
h.Normalization = 'probability';
h.BinWidth = 0.05;
h.FaceColor = 'b';
h.EdgeColor = 'none';
xlabel('ROL of insula minus vlPFC (s)')
ylabel('Probability')
xline(0,'linewidth',2)
title('all conditions')

figure('Position', [0 0 400 400])
h = histogram(ROL_asian);
h.Normalization = 'probability';
h.BinWidth = 0.05;
h.FaceColor = cdcol(1,:);
h.EdgeColor = 'none';
xlabel('ROL of insula minus vlPFC (s)')
ylabel('Probability')
xline(0,'linewidth',2)
title('asian')

figure('Position', [0 0 400 400])
h = histogram(ROL_black);
h.Normalization = 'probability';
h.BinWidth = 0.05;
h.FaceColor = cdcol(2,:);
h.EdgeColor = 'none';
xlabel('ROL of insula minus vlPFC (s)')
ylabel('Probability')
xline(0,'linewidth',2)
title('black')

figure('Position', [0 0 400 400])
h = histogram(ROL_white);
h.Normalization = 'probability';
h.BinWidth = 0.05;
h.FaceColor = cdcol(3,:);
h.EdgeColor = 'none';
xlabel('ROL of insula minus vlPFC (s)')
ylabel('Probability')
xline(0,'linewidth',2)
title('white')


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
