load('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/cell_of_44_race_cases_tables.mat');%
coords_MNI152 = [];
for ci = 1:length(T)
    coords305 = [T{ci}.IELVis_coord_1,T{ci}.IELVis_coord_2,T{ci}.IELVis_coord_3];
    