%% this is the code to creat a cell with tables

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
sbj_names_IELVis_all = {'PT020';'PT021';'PT022';'PT023';'PT024';'PT025';'PT026';'PT027';'PT028';'PT029';'PT030'...
    ;'PT031';'PT032';'PT033';'PT034';'PT035';'PT037';'PT038';'PT039';'PT040';'PT041';'PT042';'PT043';'PT044'...
    ;'PT045';'PT046';'PT047';'PT049';'PT050';'PT051';'PT052';'PT053';'PT055';'PT058';'PT060';'PT062';'S17_114'...
    ;'S17_116';'S17_118';'S19_145';'S20_148';'S20_149';'S20_150';'S20_152'};


%make a specific selection of cohort
sbj_names = sbj_names_all(1:end);sbj_names_IELVis = sbj_names_IELVis_all(1:end);%China
%sbj_names = sbj_names_all(37:end);sbj_names_IELVis = sbj_names_IELVis_all(37:end);%Stanford
fsDir='/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes';




%% Visit each excel table, add a name column, and concatenate them into a cell
T = cell(size(sbj_names,1), 1);
for i = 1:length(sbj_names)
    %cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])%plz adjust accordingly to your ecosystem
    T{i} = readtable(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i} '/' sbj_names{i} '_stats.xlsx']);
    sbj_name_channame = cell(size(T{i},1),1);
    for j = 1:size(T{i},1)
        sbj_name_channame{j} = [sbj_names{i},'-',T{i}.glv_channame{j}];
    end
    T{i}.sbj_name_channame = sbj_name_channame;
    % add MNI coordinate from Aaron's codes
    cfg = [];
    elecFname=fullfile([fsDir filesep sbj_names_IELVis{i} filesep 'elec_recon' filesep sbj_names_IELVis{i} '.electrodeNames']);
    elecInfo=csv2Cell(elecFname,' ',2);
    cfg.elecNames=append(elecInfo(:,1),elecInfo(:,3));%take care of the cfg.elecNames
    
    
    for ei = 1:length(cfg.elecNames)
        if cfg.elecNames{ei}(end) == 'R'
            cfg.isLeft(ei) = 0;
        else
            cfg.isLeft(ei) = 1;
        end
    end
    
    
    sub_coords=importdata([fsDir filesep sbj_names_IELVis{i} filesep 'elec_recon' filesep sbj_names_IELVis{i} '.PIAL']);
    
    cfg.elecCoord = sub_coords.data;
    
    [avgCoords, elecNames, isLeft]=sub2AvgBrain(sbj_names_IELVis{i},cfg);
    
    [Lia,Locb] = ismember(T{i}.LEPTO_coord_1,sub_coords.data(:,1));%the index connect pedro's and IELVis code
    T{i}.IELVis_coord_1 = zeros(size(T{i},1),1);
    T{i}.IELVis_coord_2 = zeros(size(T{i},1),1);
    T{i}.IELVis_coord_3 = zeros(size(T{i},1),1);
    for j = 1:size(T{i},1)
        if Locb(j)~=0
            T{i}.IELVis_coord_1(j) = avgCoords(Locb(j),1);
            T{i}.IELVis_coord_2(j)= avgCoords(Locb(j),2);
            T{i}.IELVis_coord_3(j)= avgCoords(Locb(j),3);
        else
            T{i}.IELVis_coord_1(j) = nan;
            T{i}.IELVis_coord_2(j) = nan;
            T{i}.IELVis_coord_3(j) = nan;
        end
    end
end
%Creat another table with rows of specific cohorts and column of specific anatomical
%structures
save('cell_of_44_race_cases_tables.mat','T')