% the following codes is for check the EDF time 
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%home

project_names = {'race_encoding_simple','race_faces'};


sbj_names = {'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
        ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58'};

% center = googleSheet.center{sbj_number};
timedata = zeros(22,1);


for i = 1:length(sbj_names)
    sbj_name = sbj_names{i};
    for j = 1:length(project_names)
        project_name = project_names{j};
        
        
        if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
            center = 'China';
        else
            center = 'Stanford';
        end
        
        block_names = BlockBySubj(sbj_name,project_name);
        dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
        [fs_iEEG, fs_Pdio, data_format] = GetFSdataFormat(sbj_name, center);
        
        
        ref_chan = [];
        epi_chan = [];
        empty_chan = [];
        
        hdr = SaveDataNihonKohden_chao(sbj_name, project_name, block_names, dirs, ref_chan, epi_chan, empty_chan);%important here downsampling in here
        
        if strcmp(project_name,'race_encoding_simple')
            t1 = hdr.starttime;
            t1(3) = ':';
            t1(6) = ':';
            t11 = datevec(datenum(t1));
        else
            t2 = hdr.starttime;
            t2(3) = ':';
            t2(6) = ':';
            t22 = datevec(datenum(t2));
            
            indx = find(ismember(sbj_names,sbj_name));
            timedata(indx) = etime(t22,t11);
        end
        
        
    end
end


load cdcol.mat
clear g
g=gramm('x',timedata(:,1)./60./60,'y',timedata(:,2),'color',task(:));
g.geom_point();
g.set_names('x','relative time(h)','y','cases','color','task');
g.set_title('the relative time of race_encoding and race_faces');
g.set_color_options('map',[cdcol.blue_jeans;cdcol.rose_pink]);
g.axe_property('Ylim',[0 22]);
g.set_point_options('base_size',8);
g.set_text_options('base_size',22,...
    'label_scaling',1,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1);
g.draw();



