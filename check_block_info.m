% the following codes is for check the block info for a specificic task
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_personal-master/'))%change here
addpath(genpath('/Users/chao/Documents/Stanford/code/lbcn_preproc-master/'))%change here
[server_root, comp_root, code_root] = AddPaths('Chao_iMAC');%change here

project_name = 'Grad_CPT';%change here


sbj_names = {'C18_30'};%change here


block_data = cell(length(sbj_names),2);


for i = 1:length(sbj_names)
    sbj_name = sbj_names{i};
    
    
    
    if contains(sbj_name,'C17')||contains(sbj_name,'C18')||contains(sbj_name,'C19')
        center = 'China';
    else
        center = 'Stanford';
    end
    
    block_names = BlockBySubj(sbj_name,project_name);
    
    block_data{i,1} = length(block_names);
    
    block_times = zeros(1,length(block_names));
    
    dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);
    [fs_iEEG, fs_Pdio, data_format] = GetFSdataFormat(sbj_name, center);
    
    for j = 1:length(block_names)
        ref_chan = [];
        epi_chan = [];
        empty_chan = [];
        bn = block_names{j}
        fn = sprintf('%s/originalData/%s/global_%s_%s_%s.mat',dirs.data_root,sbj_name,project_name,sbj_name,bn);
        load(fn,'globalVar');
        
        data_dir = [dirs.original_data '/' sbj_name '/' bn]; 
        fname =  [dirs.original_data '/' sbj_name '/' bn '/' bn '.edf'];
        
        if strcmp(globalVar.center, 'Stanford')
            [hdr, D] = edfread(fname);
        elseif strcmp(globalVar.center, 'China')
            if strcmp(sbj_name, 'C18_33')
                [hdr,D] = edfread_China(fname,'targetsignals',[1:105]);
            elseif  strcmp(sbj_name, 'C18_35')
                [hdr,D] = edfread_China(fname,'targetsignals',[1:186]);
            elseif  strcmp(sbj_name, 'C18_38')
                [hdr,D] = edfread_China(fname,'targetsignals',[1:146]);
            elseif  strcmp(sbj_name, 'C18_42')
                [hdr,D] = edfread_China(fname,'targetsignals',[1:123]);
            elseif  strcmp(sbj_name, 'C18_32')
                [hdr,D] = edfread_China(fname,'targetsignals',[1:163]);
            else
                [hdr, D] = edfread_China(fname);
            end
        else
        end
        
        block_times(1,j) = hdr.records * hdr.duration/60;
         % hdr.records = number of chuncks 
         % hdr.duration = duration of each chunck 
    end
    
    block_data{i,2} = block_times;
end

block_numbers = vertcat(block_data{:,1});

block_times = vertcat(block_data{:,2});
block_times = block_time(:);

fprintf('subjects performed between %d-%d runs (mean %d) (%d-%d minutes, mean %d minutes ) of %s ',...
    min(block_numbers),max(block_numbers),mean(block_numbers),min(block_times),max(block_times),mean(block_times),project_name);

