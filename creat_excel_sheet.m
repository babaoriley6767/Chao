    %% write excel sheet based on subjvar
    % this part is about build excel files
    sbj_names = {'C17_20';'C17_21';'C18_22';'C18_23';'C18_24';'C18_25';'C18_26';'C18_27';'C18_28';'C18_29';'C18_30'...
        ;'C18_31';'C18_32';'C18_33';'C18_34';'C18_35';'C18_37';'C18_38';'C18_39';'C18_40';'C18_41';'C18_42';'C18_43';'C18_44'...
        ;'C18_45';'C18_46';'C18_47';'C18_49';'C19_50';'C19_51';'C19_52';'C19_53';'C19_55';'C19_58';'C19_60';'C19_62';'S17_114_EB'...
        ;'S17_116_AA';'S17_118_TW';'S20_148_SM';'S20_149_DR';'S20_150_CM';'S20_152_HT';'S19_145_PC'};
    for i = 1:length(sbj_names)
        cd(['/Volumes/CHAO_IRON_M/data/neuralData/originalData/' sbj_names{i}])
        load(['subjVar_' sbj_names{i} '.mat'])
        elinfo_link = subjVar.elinfo;
        [file path] = uigetfile;%select the GlobalVar
        load([path file]);
        glv_channame = globalVar.channame';
        if size(elinfo_link,1) == size(glv_channame,1)
            disp('Glv and freesurfer match perfect')
            glv_index = [1:globalVar.nchan]';
        elseif size(elinfo_link,1) > size(glv_channame,1)
            disp('Glv naming is shorter than freesurfer!!! Adjust the .XLSX')
            for j = size(glv_channame,1)+1:size(elinfo_link,1)
                glv_channame{j} = 'empty';
            end
            glv_index = [1:size(elinfo_link,1)]';
        else
            warning('the lenght of GlobalVar is larger than the Freesurfer')
            if strcmp(sbj_names{i},'C17_21')
                chan_in_glv = [1:19,21:131];
                glv_channame = glv_channame(chan_in_glv);
                glv_index = [1:globalVar.nchan]';
                glv_index = glv_index(chan_in_glv);
            elseif strcmp(sbj_names{i},'C18_29')
                chan_in_glv = [1:19,21:159];
                glv_channame = glv_channame(chan_in_glv);
                glv_index = [1:globalVar.nchan]';
                glv_index = glv_index(chan_in_glv);
            elseif strcmp(sbj_names{i},'C18_31')
                chan_in_glv = [1:126];
                glv_channame = glv_channame(chan_in_glv);
                glv_index = [1:globalVar.nchan]';
                glv_index = glv_index(chan_in_glv);
            elseif strcmp(sbj_names{i},'C18_45')
                chan_in_glv = [1:236];
                glv_channame = glv_channame(chan_in_glv);
                glv_index = [1:globalVar.nchan]';
                glv_index = glv_index(chan_in_glv);
            elseif strcmp(sbj_names{i},'C19_62')
                chan_in_glv = [1:128];
                glv_channame = glv_channame(chan_in_glv);
                glv_index = [1:globalVar.nchan]';
                glv_index = glv_index(chan_in_glv);
            else
                
            end
        end
        elinfo_link = addvars(elinfo_link,glv_channame,'Before','FS_vol');
        elinfo_link = addvars(elinfo_link,glv_index,'Before','FS_vol');
        writetable(elinfo_link, [subjVar.sbj_name '_new.xlsx'] );
        macopen([subjVar.sbj_name '_new.xlsx'])
        prompt = 'Do you align the FS_label and glv_channame,if yes click y/';
        ID = input(prompt,'s');
        if strcmp(ID, 'y')
            disp('good to go');
        else
            warning(['please check the subjvar' sbj_names{i}])
            RETURN;
        end
    end
    