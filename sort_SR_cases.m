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
anat = {'INSULA'};anat_name = 'INSULA';
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
anat = {'SFS','SFG'};anat_name = 'SFS SFG';
anat = {'ACC','MCC'};anat_name = 'ACC MCC';


blocks = {'race','race_faces','VTC','gradCPT'};
side = 'none';%'L','R','none'
channame_4 = cell(1,4);

for project = 1:length(blocks)
    if strcmp(blocks{project},'race')
        indxcohort = 1:36;sbj_names = sbj_names_all(indxcohort);%China
        [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');
        channame_4{project} = channame;
    elseif strcmp(blocks{project},'race_faces')
        sbj_names_all_race_faces = {'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
            ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_52';'C19_53';'C19_55';'C19_58'};
        indxcohort = ismember(sbj_names_all,sbj_names_all_race_faces);sbj_names = sbj_names_all(indxcohort);%China
        [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');
        channame_4{project} = channame;
    elseif strcmp(blocks{project},'VTC')
        sbj_names_all_VTC = {'C17_20';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
            ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
            ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_55';'C19_60';'C19_62'};       
        indxcohort = ismember(sbj_names_all,sbj_names_all_VTC);sbj_names = sbj_names_all(indxcohort);%China
        [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');
        channame_4{project} = channame;
    else
        
        sbj_names_all_gradCPT = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_29';'C18_30'...
            ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
            ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62'};
        indxcohort = ismember(sbj_names_all,sbj_names_all_gradCPT);sbj_names = sbj_names_all(indxcohort);%China
        [channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');
        channame_4{project} = channame;
    end
end
intersect(intersect(intersect(channame_4{1},channame_4{2}),channame_4{3}),channame_4{4})


%% Visit each excel table, add a name column, and concatenate them into a cell
