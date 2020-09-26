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
% sbj_names = sbj_names_all(1:36);%China
% sbj_names = sbj_names_all(37:end);%Stanford
sbj_names = sbj_names_all;%all



% define the the abbreviations of kinds of brian structures
anat_all = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','WM','OUT','EC',...
    'FOP','POP','TOP','EMPTY','PARL','LESION','INSULA','BASAL'};
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);


%make a specific selection of anatomical structures
anat = {'INSULA','ACC','MCC','FG','IFS','IFG','SFS','SFG','AMY'};



%%  extract the table information through subjects and anat
sbj_name = [];
sbj_name_insheet = [];
FS_label = [];
label = [];
LvsR = [];

T = cell(size(sbj_names,1), 1);
for i = 1:length(sbj_names)
    T{i} = readtable(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i} '/' sbj_names{i} '_stats.xlsx']);
    sbj_name_insheet = cell(size(T{i},1),1);
    for j = 1:size(T{i},1)
        sbj_name_insheet{j} = [sbj_names{i}];
    end
    T{i}.sbj_name = sbj_name_insheet;
end
for i = 1:length(sbj_names)
    for j = 1:length(anat)

         
        idx1 = strcmp(T{i}.label,anat{j});
        idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
        idx = idx1 & idx2;
        
        FS_label = [FS_label;T{i}.FS_label(idx)];
        label = [label;T{i}.label(idx)];
        LvsR = [LvsR;T{i}.LvsR(idx)];
        sbj_name = [sbj_name;T{i}.sbj_name(idx)];
    end
end
%% build a table and make a sheet
T_any_activation = table(sbj_name,FS_label,LvsR,label);

writetable(T_any_activation, 'any_activation_and_anat.xlsx' );

    %% move the excel files to a folder.
    sbj_names = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
        ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
        ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
        ;'S17_116_AA';'S17_118_TW';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};
    for i = 1:length(sbj_names)
        cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])
        copyfile ([sbj_names{i} '_stats.xlsx'] ,'/Users/chao/Desktop/Clara_localization')
    end
