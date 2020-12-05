function [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,stats,side)
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
end