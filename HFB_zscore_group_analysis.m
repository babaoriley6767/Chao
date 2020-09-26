clear;clc;
cd /Users/chao/Desktop/Project_in_Stanford/RACE/4_working_data/HFB_zscore

pathformation = dir;
pathnumber = numel(pathformation);

group_index = cumsum([0 2 1 2 3 5]);% the index for the final


for ni = 3:pathnumber
    load(pathformation(ni).name)
    data = hfb_mean_data;
    switch data.sbj_name
        case 'C18_23'
            index = {'X1','X2'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT023-',index{i},'L'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
        case 'C18_28'
            index = {'Q2'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT028-',index{i},'R'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
        case 'C18_29'
            index = {'X''6','X''7'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT029-',index{i}(1),index{i}(3),'L'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
        case 'C18_37'
            index = {'X4','X5','X6'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT037-',index{i},'R'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
        case 'C18_49'
            index = {'X4','X5','X6','X7','X8'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT049-',index{i},'R'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
        case 'C19_62'
            index = {'R1','X6'};
            for i = 1:length(index)
                grouphfb.name{i+group_index(ni-2),1} = ['PT062-',index{i},'R'];
                index_of_contact = find(strcmp(data.label,index{i}));
                grouphfb.hfb{i+group_index(ni-2),1} = data.meandata(index_of_contact,:);
            end
    end
end


save('/Volumes/CHAO_IRON_M/iELVis_working/Plot_Elelctrodes/grouphfb.mat','grouphfb')