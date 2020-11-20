
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
anat = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','EC',...
    'FOP','POP','TOP','PARL','INSULA','BASAL'};anat_name='all available';
anat = {'WM','OUT','EMPTY','LESION'};;anat_name='exclude';

side = 'none';%'L','R','none'



anat = {'STS'};anat_name='STS';
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);
%% define the loop
anat = cell(9,1);
anat_name = cell(9,1);
anat{1} = {'INSULA'};anat_name{1} = 'INSULA';
anat{2} = {'IFS','IFG'};anat_name{2} = 'IFS IFG';%
anat{3} = {'SFS','SFG'};anat_name{3} = 'SFS SFG';
anat{4} = {'ACC','MCC'};anat_name{4} = 'ACC MCC';
anat{5} = {'AMY'};anat_name{5} = 'amygdala';
anat{6} = {'HIPPO A'};anat_name{6} = 'HIPPO A';
anat{7} = {'HIPPO M'};anat_name{7} = 'HIPPO M';
anat{8} = {'HIPPO P'};anat_name{8} = 'HIPPO P';
anat{9} = {'PHG'};anat_name{9} = 'PHG';
%% Visit each excel table, add a name column, and concatenate them into a cell
for k = 1:length(anat)
    load('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/cell_of_44_race_cases_tables.mat');%if there is any change of the excel sheet,
    %then this need to update,go to 'Creat_cell_of_tables.mat'
    T = T(indxcohort,1);
    
    %Creat another table with rows of specific cohorts and column of specific anatomical
    %structures
    sz = [size(sbj_names,1) size(anat{k},2)];
    varTypes = cell(1,size(anat{k},2));
    varTypes(:) = {'cell'};
    T2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',anat{k},'RowNames',sbj_names);
    %put the glv_index into each space of the table as a vector
    if isempty(side)||strcmp(side,'none')
        for i = 1:length(sbj_names)
            for j = 1:length(anat{k})
                idx1 = strcmp(T{i}.label,anat{k}{j});
                idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation/T{i}.group_diff %default
                idx = idx1 & idx2;
                T2{sbj_names{i},anat{k}{j}} = {T{i}.glv_index(idx)'};
            end
        end
    else
        for i = 1:length(sbj_names)
            for j = 1:length(anat{k})
                idx1 = strcmp(T{i}.label,anat{k}{j});
                idx2 = T{i}.any_activation;%or idx2 = T{i}.any_activation/T{i}.all_trials_activation
                idx3 = strcmp(T{i}.LvsR,side);
                idx = idx1 & idx2 & idx3;
                T2{sbj_names{i},anat{k}{j}} = {T{i}.glv_index(idx)'};
            end
        end
    end
    %Since there may be empty electrodes in the stats sheet, Glv_index may be str. Here, all Glv_index in T2 will be transformed into vetor
    for i = 1:length(sbj_names)
        for j = 1:length(anat{k})
            if iscell(T2{sbj_names{i},anat{k}{j}}{:})
                T2{sbj_names{i},anat{k}{j}}{:} = str2double(T2{sbj_names{i},anat{k}{j}}{:});
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
    
behv = readtable(['/Users/chao/Desktop/Project_in_Stanford/01_RACE/4_working_data/Behavioral_accuracy/results_summary.xlsx']);
behv_indx = ismember(behv.Chao_patient_ID_in_server,T3.Properties.RowNames);
CatAcc_SR = behv.Race_CatAcc_SR(behv_indx);
CatAcc_OR = behv.Race_CatAcc_OR(behv_indx);
% disp(['the mean of SR accuracy in ' anat{k}{:} ' is ' num2str(round(mean(CatAcc_SR)*1000)/10) ' and the std is ' num2str(round(std(CatAcc_SR)*1000)/10)]);
% disp(['the mean of OR accuracy in ' anat{k}{:} ' is ' num2str(round(mean(CatAcc_OR)*1000)/10) ' and the std is ' num2str(round(std(CatAcc_OR)*1000)/10)]);
disp(['the mean of total accuracy in ' anat{k}{:} ' is ' num2str(round(mean((CatAcc_SR+CatAcc_OR)./2)*1000)/10) ' and the std is ' num2str(round(std((CatAcc_SR+CatAcc_OR)./2)*1000)/10)]);
end





