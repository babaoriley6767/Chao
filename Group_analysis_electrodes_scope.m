%get ready with the addpath,plz adjust accordingly
% this is for the question of JNS " Effects in the Insula were probed more deeply than in other structures, including examining brain-behavior correlations and anatomical distribution, but no justification was given for concentrating on the insula in this way."
addpath(genpath('/Users/tony/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/tony/Documents/Stanford/code/lbcn_preproc-master/'))
%addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master'))% this is a matlab based graph toolbox which is similar to the ggplot2
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

% define the cohort
sbj_names_all = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
    ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
    ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
    ;'S17_116_AA';'S17_118_TW';'S19_145_PC';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};

%make a specific selection of cohort
indxcohort = 1:36;%China


sbj_names = sbj_names_all(indxcohort);%China

anat_all = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','EC',...
    'FOP','POP','TOP','PARL','INSULA','BASAL'};




side = 'none';%'L','R','none'

%%

for q = 1:length(anat_all)
    
anat = {anat_all{1,q}};
anat_name = anat_all{1,q};



%% Visit each excel table, add a name column, and concatenate them into a cell
% [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');  
stats = 'group_diff';
load('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/cell_of_44_race_cases_tables.mat');%if there is any change of the excel sheet, 
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

len = cellfun('prodofsize', T3.anat);
n = sum(len);
%%
 disp([anat ' has ' num2str(height(T3)) ' patiens and ' num2str(n) ' group_diff electrodes'] )
 
 
 
 
 end