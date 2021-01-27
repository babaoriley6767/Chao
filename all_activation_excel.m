clear% %
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home
project_name = 'race_encoding_simple';
%% the following part is about index in the excel (activated or not)
sbj_names = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
        ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
        ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
        ;'S17_116_AA';'S17_118_TW';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT'};
    for i = 1:length(sbj_names)
        if contains(sbj_names{i},'C17')||contains(sbj_names{i},'C18')||contains(sbj_names{i},'C19')
            center = 'China'; 
        else
            center = 'Stanford';
        end
        dirs = InitializeDirs(project_name, sbj_names{i}, comp_root, server_root, code_root);
        block_names = BlockBySubj(sbj_names{i},project_name);
        cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])
        load([dirs.data_root,filesep,'OriginalData',filesep,sbj_names{i},filesep,'global_',project_name,'_',sbj_names{i},'_',block_names{1},'.mat'])
        T = readtable([sbj_names{i},'_stats.xlsx']);
        
        stats_params = genStatsParams(project_name);
        stats_params.paired = true; 
        stats_params.all_trials = true; condNames = 'all_trials';conds = {'all_trials'};
       
        [p,p_fdr,sig,greater] = permutationStatsAll(sbj_names{i},project_name,block_names,dirs,[],'stimlock_bl_corr',condNames,conds,'Band','HFB',stats_params);
        sig_AT_bl = intersect(find(sig==1),find(greater==1));
        sig_non_AT_bl = find(sig==0);
        
        
        if ~iscell(T.glv_index)
            all_trials_activation = ismember(T.glv_index,sig_AT_bl);
        else
            sig_AT_bl_str = cell(length(sig_AT_bl),1);
            for j = 1:length(sig_AT_bl)
                sig_AT_bl_str{j} = num2str(sig_AT_bl(j),'%2d');
            end
              all_trials_activation = ismember(T.glv_index,sig_AT_bl_str);
            %any_activation = ismember(T.glv_index,mat2cell(num2str(ind_abw,'%2d'),ones(1,length(ind_abw))));
        end
        T.all_trials_activation = all_trials_activation;
        writetable(T, [sbj_names{i} '_stats.xlsx'] );
    end
