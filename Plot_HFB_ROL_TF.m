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
indxcohort = 1:36;%China
% indxcohort = [37:44];%Stanfordc
% indxcohort = [1:44];%2 Centers

sbj_names = sbj_names_all(indxcohort);%China


% define the the abbreviations of kinds of brian structures
anat_all = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','WM','OUT','EC',...
    'FOP','POP','TOP','EMPTY','PARL','LESION','INSULA','BASAL'};

%make a specific selection of anatomical structures
anat = {'HIPPO A','HIPPO M','HIPPO P','PHG'};anat_name = 'MTLE';
anat = {'HIPPO A'};anat_name = 'HIPPO A';
anat = {'HIPPO M'};anat_name = 'HIPPO M';
anat = {'HIPPO P'};anat_name = 'HIPPO P';
anat = {'PHG'};anat_name = 'PHG';
anat = {'INSULA'};anat_name = 'INSULA';
anat = {'ACC','MCC'};anat_name = 'ACC MCC';
anat = {'ACC'};anat_name = 'ACC';
anat = {'MCC'};anat_name = 'MCC';
anat = {'FG'};anat_name = 'FG';
anat = {'FG','OTS','CS'};anat_name = 'FG OTS CS';
anat = {'IFS','IFG'};anat_name = 'IFS IFG';%
anat = {'SFS','SFG'};anat_name = 'SFS SFG';
anat = {'OFC'};anat_name = 'OFC';
anat = {'AMY'};anat_name = 'amygdala';
anat = {'SMA','PARL'};anat_name = 'SMA and PARL';
anat = {'SOG','MOG','IOG'}; anat_name= 'occipital';
anat = {'FOP','POP','TOP'}; anat_name= 'operculum';
anat = {'PCG','CG','POF','CF','PCC'}; anat_name='posterior medial';
anat = {'SPL','IPL','IPS'};anat_name='parieltal area';
anat = {'POSTCENTRAL G','PRECENTRAL G'};anat_name='central area';
anat = {'POSTCENTRAL G'};anat_name='central area';
anat = {'STG','STS','MTG','ITS','ITG'};anat_name='temporal area lateral';
anat = {'PCG'};anat_name='PCG';
anat = {'MTG'};anat_name = 'MTG';
anat = {'IFS','IFG',}
anat = {'INSULA','IFS','IFG','SFS','SFG','ACC','MCC'};anat_name = 'ROI';
anat = {'SFG','SFS','MFG','IFS','IFG','OFC','MPG','SMA','VMPFC','ACC','MCC','PCC','STG','STS','MTG','ITS','ITG','AMY','HIPPO A','HIPPO M','HIPPO P'...
    ,'TP','OTS','FG','CS','PHG','PRECENTRAL G','POSTCENTRAL G','SPL','IPL','IPS','PCG','CG','POF','CF','LG','SOG','MOG','IOG','EC',...
    'FOP','POP','TOP','PARL','INSULA','BASAL'};anat_name='all available';
anat = {'WM','OUT','EMPTY','LESION'};;anat_name='exclude';

side = 'none';%'L','R','none'
site_pick = 'all sites per subj';%all sites per subj

site_pick = 'One site  per subj';%One touch per person


anat = {'STS'};anat_name='STS';
% Check on the meaning of the abbreviations
anat_displ =  importdata('/Users/chao/Documents/Stanford/code/lbcn_personal-master/Chao/anat_abbreviation.txt');%pls select a directory to store the 
disp(anat_displ);


%% choose sites from sheet
% Visit each excel table, add a name column, and concatenate them into a cell
[channame,T3] = group_analysis_part1(sbj_names,indxcohort,anat,'group_diff','none');   
% paticular sites
switch anat_name
    case 'INSULA'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,2,5,6],:);
            T3.anat{3} = T3.anat{3}(3);
        else
        end
    case 'IFS IFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([5,10],:);
            T3.anat{1} = T3.anat{1}(5);
        else
        end
    case 'SFS SFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,4,6],:);
            T3.anat{2} = T3.anat{2}(3);
            T3.anat{3} = T3.anat{3}(1);
        else
        end
end
disp(channame);

%load sheet and vertcat sheet
load('/Users/tony/Documents/Stanford/code/lbcn_personal-master/Chao/mat_files/group_diff_encoding_HFB_reall_behav_sheet.mat')
channame_sheet = cell(length(HFB_behav),1);
for i = 1:length(HFB_behav)
    channame_sheet{i} = HFB_behav{i}.channame{1};
end
HFB_behav = HFB_behav(1,ismember(channame_sheet,channame));
HFB_behav_stats = vertcat(HFB_behav{:});

%% 
%remove outlier of ROL_onset (modified z-score)

zscorethresh = 3;

idx_clean = HFB_behav_stats.clean_trial==1;


idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);


idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);


ROL_asian = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_asian);
ROL_black = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_black);
ROL_white = HFB_behav_stats.ROL_onsets(idx_ROL_nan&idx_white);


mZ_column = nan(height(HFB_behav_stats),1);%% zscore
for i = 1:height(HFB_behav_stats)
    if HFB_behav_stats.test_race(i) == 1 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_asian)) ./ mad(ROL_asian,1);
    elseif  HFB_behav_stats.test_race(i) == 2 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_black)) ./ mad(ROL_black,1);
    elseif  HFB_behav_stats.test_race(i) == 3 &&  ~isnan(HFB_behav_stats.ROL_onsets(i))
        mZ_column(i,1) = norminv(.75)*(HFB_behav_stats.ROL_onsets(i)-median(ROL_white)) ./ mad(ROL_white,1);
    else
    end
end

mZ_indx = ~(mZ_column>zscorethresh);% key parameter
ROL_asian_mZ = mZ_column(idx_ROL_nan&idx_asian&mZ_indx);% After zscore outlier exclusion 
ROL_black_mZ = mZ_column(idx_ROL_nan&idx_black&mZ_indx);
ROL_white_mZ = mZ_column(idx_ROL_nan&idx_white&mZ_indx);


%%
% plot the figure
%Asian
figure(1), clf
subplot(311), hold on
plot(ROL_asian,'k^','markerfacecolor','w','markersize',12);
set(gca,'xtick',[]), box off
ylabel('ROL')
title('mZscore Asian')

subplot(312), hold on
plot(mZ_column(idx_ROL_nan&idx_asian),'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore outlier')

subplot(313), hold on
plot(ROL_asian_mZ,'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore')

%Black
figure(2), clf
subplot(311), hold on
plot(ROL_black,'k^','markerfacecolor','w','markersize',12);
set(gca,'xtick',[]), box off
ylabel('ROL')
title('mZscore Black')

subplot(312), hold on
plot(mZ_column(idx_ROL_nan&idx_black),'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore outlier')

subplot(313), hold on
plot(ROL_black_mZ,'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore')

%White
figure(3), clf
subplot(311), hold on
plot(ROL_white,'k^','markerfacecolor','w','markersize',12);
set(gca,'xtick',[]), box off
ylabel('ROL')
title('mZscore White')

subplot(312), hold on
plot(mZ_column(idx_ROL_nan&idx_white),'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore outlier')

subplot(313), hold on
plot(ROL_white_mZ,'k^','markerfacecolor','w','markersize',12);
plot(get(gca,'xlim'),[1 1]*zscorethresh,'--','color','r')
set(gca,'xtick',[]), box off
ylabel('Zscore')

%%
%time-trial-hfb plot of ROL onset asian black and white

HFB_range = dsearchn(data_time',[0 1]');

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian&mZ_indx,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black&mZ_indx,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white&mZ_indx,HFB_range(1):HFB_range(2));

data_asian_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian&mZ_indx,:);
[data_asian_ROL_onset_s,asian_sort] = sort(data_asian_ROL_onset);
data_asian_s = data_asian(asian_sort,:);
data_black_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black&mZ_indx,:);
[data_black_ROL_onset_s,black_sort] = sort(data_black_ROL_onset);
data_black_s = data_black(black_sort,:);
data_white_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white&mZ_indx,:);
[data_white_ROL_onset_s,white_sort] = sort(data_white_ROL_onset);
data_white_s = data_white(white_sort,:);

%time-trial-hfb plot of ROL onset asian black and white(without modified Zscore outlier exclusion) 
HFB_range = dsearchn(data_time',[0 1]');

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_onset_s,asian_sort] = sort(data_asian_ROL_onset);
data_asian_s = data_asian(asian_sort,:);
data_black_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_onset_s,black_sort] = sort(data_black_ROL_onset);
data_black_s = data_black(black_sort,:);
data_white_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_onset_s,white_sort] = sort(data_white_ROL_onset);
data_white_s = data_white(white_sort,:);

%% plot time trial hfb results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_onset_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
% title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Asian'])
title('Asian')
set(gca,'fontsize',28)
ylabel('trials')
h = colorbar;
ylabel(h, 'Z-scored power')

subplot(3,1,2)
contourf(data_time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_onset_s,size(data_black,1):-1:1,'k--','LineWidth',3)
% title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Black'])
title('Black')
set(gca,'fontsize',28)
ylabel('trials')
h = colorbar;
ylabel(h, 'Z-scored power')

subplot(3,1,3)
contourf(data_time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_onset_s,size(data_white,1):-1:1,'k--','LineWidth',3)
% title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, White'])
title('White')
set(gca,'fontsize',28)
ylabel('trials')
h = colorbar;
ylabel(h, 'Z-scored power')
xlabel('Time(S)')


data_stats = [data_asian_ROL_onset;data_black_ROL_onset;data_white_ROL_onset];
data_groups = [ones(length(data_asian_ROL_onset),1);ones(length(data_black_ROL_onset),1)*2;ones(length(data_white_ROL_onset),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups)
multcompare(stats,'ctype','lsd')% the anova stats of ROL
disp(['asian value ' num2str(mean(data_asian_ROL_onset)*1000) ' ± ' num2str(std(data_asian_ROL_onset)*1000)])
disp(['black value ' num2str(mean(data_black_ROL_onset)*1000) ' ± ' num2str(std(data_black_ROL_onset)*1000)])
disp(['white value ' num2str(mean(data_white_ROL_onset)*1000) ' ± ' num2str(std(data_white_ROL_onset)*1000)])

%%
asian_hfb = HFB_behav_stats.("HFB 0.10~1.00")(idx_clean&idx_ROL_nan&idx_asian&mZ_indx,:);
black_hfb = HFB_behav_stats.("HFB 0.10~1.00")(idx_clean&idx_ROL_nan&idx_black&mZ_indx,:);
white_hfb = HFB_behav_stats.("HFB 0.10~1.00")(idx_clean&idx_ROL_nan&idx_white&mZ_indx,:);



rho = cell(1,3);
pval = cell(1,3);
[rho{1},pval{1}] = corr(data_asian_ROL_onset,asian_hfb);
[rho{2},pval{2}] = corr(data_black_ROL_onset,black_hfb);
[rho{3},pval{3}] = corr(data_white_ROL_onset,white_hfb);



%%

ROL_onset_column = data_stats;
HFB_column = [asian_hfb;black_hfb;white_hfb];
group_column = data_groups;

group_column_cell = cell(length(group_column),1);
group_column_cell(group_column==1,1) = deal({'Asian'});
group_column_cell(group_column==2,1) = deal({'Black'});
group_column_cell(group_column==3,1) = deal({'White'});

clear g
close all
figure('Position', [1500 500 700 700]),clf
g(1,1)=gramm('x',ROL_onset_column,'color',group_column_cell);
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','');
% g(1,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',45); %histogram
g(1,1).stat_bin('normalization','cumcount','geom','line')
g(1,1).axe_property('XTickLabel','','Fontsize',15); % We deactivate tht ticks


%Create a scatter plot
g(2,1)=gramm('x',ROL_onset_column,'y',HFB_column,'color',group_column_cell);
g(2,1).set_names('x','ROL onset time','y','HFB zscore','color','Group');
g(2,1).geom_point(); %Scatter plot
g(2,1).stat_glm();
g(2,1).set_point_options('base_size',6);
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend_pos',[0.83 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(2,1).axe_property('Ygrid','on','Fontsize',15); 
g(2,1).legend_axe_handle.FontSize = 25;
% set(g(2,1).legend_axe_handle.FontSize,30);
% h = legend('x','y','z');
%     set(h,'FontSize',28);


%Create y data histogram on the right
g(3,1)=gramm('x',HFB_column,'color',group_column_cell);
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','');
g(3,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',45); %histogram
% g(3,1).stat_bin('normalization','cumcount','geom','stairs')
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');
g(3,1).axe_property('Fontsize',15); 

g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('1s HFB correlate ROL','Fontsize',15);
g.set_color_options('map',plot_params.col);
% g.set_color_options('map','d3_10');
g.draw();
% set(g(2,1).legend_axe_handle.FontSize,30);

%% violin plot of ROL
ROL_onset_column = data_stats;
% HFB_column = [asian_hfb;black_hfb;white_hfb];
group_column = data_groups;

group_column_cell = cell(length(group_column),1);
group_column_cell(group_column==1,1) = deal({'Asian'});
group_column_cell(group_column==2,1) = deal({'Black'});
group_column_cell(group_column==3,1) = deal({'White'});


% 
% data_stats = [ROL_OR;ROL_SR];
% group_OR = cell(length(ROL_OR),1);
% group_OR(:) = deal({'OR'});
% group_SR = cell(length(ROL_SR),1);
% group_SR(:) = deal({' SR'});%add space this data column will come to the first
% data_groups = vertcat(group_OR,group_SR);
% clear g
% 
% g(1,1)=gramm('x',group_column_cell,'y',ROL_onset_column,'color',group_column_cell,'linestyle','--');
% g(1,1).stat_violin('fill','transparent');
% g(1,1).set_title('stat_violin()');

g=gramm('x',group_column_cell,'y',ROL_onset_column,'color',group_column_cell);
g.set_names('x',[],'y','ROL(S)','color','Origin');
g.stat_violin('normalization','area','dodge',0,'fill','edge');
g.stat_boxplot('width',0.2);
g.set_color_options('map',plot_params.col);
g.axe_property('FontSize',28)


figure('Position',[100 100 600 550]);
g.draw();



%%
load example_data;
load fisheriris.mat

clear g
figure('Position',[100 100 550 550]);

%Create x data histogram on top
g(1,1)=gramm('x',meas(:,2),'color',species);
g(1,1).set_layout_options('Position',[0 0.8 0.8 0.2],... %Set the position in the figure (as in standard 'Position' axe property)
    'legend',false,... % No need to display legend for side histograms
    'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    'margin_width',[0.1 0.02],...
    'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
g(1,1).set_names('x','');
g(1,1).stat_bin('geom','overlaid_bar','fill','transparent','nbins',15); %histogram
g(1,1).axe_property('XTickLabel',''); % We deactivate tht ticks

%Create a scatter plot
g(2,1)=gramm('x',meas(:,2),'y',meas(:,1),'color',species);
g(2,1).set_names('x','Sepal Width','y','Sepal Length','color','Species');
g(2,1).geom_point(); %Scatter plot
g(2,1).set_point_options('base_size',6);
g(2,1).set_layout_options('Position',[0 0 0.8 0.8],...
    'legend_pos',[0.83 0.75 0.2 0.2],... %We detach the legend from the plot and move it to the top right
    'margin_height',[0.1 0.02],...
    'margin_width',[0.1 0.02],...
    'redraw',false);
g(2,1).axe_property('Ygrid','on'); 

%Create y data histogram on the right
g(3,1)=gramm('x',meas(:,1),'color',species);
g(3,1).set_layout_options('Position',[0.8 0 0.2 0.8],...
    'legend',false,...
    'margin_height',[0.1 0.02],...
    'margin_width',[0.02 0.05],...
    'redraw',false);
g(3,1).set_names('x','');
g(3,1).stat_bin('geom','stacked_bar','fill','all','nbins',15); %histogram
g(3,1).coord_flip();
g(3,1).axe_property('XTickLabel','');

%Set global axe properties
g.axe_property('TickDir','out','XGrid','on','GridColor',[0.5 0.5 0.5]);
g.set_title('Fisher Iris, custom layout');
g.set_color_options('map','d3_10');
g.draw();











%%
% time-trial-hfb plot of ROL peaks asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;



idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_peaks);

idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_peaks_s,asian_sort] = sort(data_asian_ROL_peaks);
data_asian_s = data_asian(asian_sort,:);

data_black_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_peaks_s,black_sort] = sort(data_black_ROL_peaks);
data_black_s = data_black(black_sort,:);

data_white_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_peaks_s,white_sort] = sort(data_white_ROL_peaks);
data_white_s = data_white(white_sort,:);

% plot results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_peaks_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_peaks_s,size(data_black,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_peaks_s,size(data_white,1):-1:1,'k--','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, White'])

data_stats = [data_asian_ROL_peaks;data_black_ROL_peaks;data_white_ROL_peaks];
data_groups = [ones(length(data_asian_ROL_peaks),1);ones(length(data_black_ROL_peaks),1)*2;ones(length(data_white_ROL_peaks),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')
%%
% time-trial-hfb plot of ROL mixed onset and peaks, asian black and white
%time-trial-hfb plot of ROL onset asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;


idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);


idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_onset_s,asian_sort] = sort(data_asian_ROL_onset);
data_asian_s = data_asian(asian_sort,:);
asian_onsets_peaks = HFB_behav_stats.ROL_peaks(asian_sort,:);

data_black_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_onset_s,black_sort] = sort(data_black_ROL_onset);
data_black_s = data_black(black_sort,:);
black_onsets_peaks = HFB_behav_stats.ROL_peaks(black_sort,:);

data_white_ROL_onset = HFB_behav_stats.ROL_onsets(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_onset_s,white_sort] = sort(data_white_ROL_onset);
data_white_s = data_white(white_sort,:);
white_onsets_peaks = HFB_behav_stats.ROL_peaks(white_sort,:);


figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_onset_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
plot(asian_onsets_peaks,size(data_asian,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_onset_s,size(data_black,1):-1:1,'k--','LineWidth',3)
plot(black_onsets_peaks,size(data_black,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_onset_s,size(data_white,1):-1:1,'k--','LineWidth',3)
plot(white_onsets_peaks,size(data_white,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL onset, White'])

data_stats = [data_asian_ROL_onset;data_black_ROL_onset;data_white_ROL_onset];
data_groups = [ones(length(data_asian_ROL_onset),1);ones(length(data_black_ROL_onset),1)*2;ones(length(data_white_ROL_onset),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')

%%
% time-trial-hfb plot of ROL mixed peaks and onsets, asian black and white
HFB_behav_stats = vertcat(HFB_behav{:});

HFB_range = dsearchn(data_all.time',[0 1]');

idx_clean = HFB_behav_stats.clean_trial==1;



idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_onsets)&~isnan(HFB_behav_stats.ROL_peaks);
%idx_ROL_nan = ~isnan(HFB_behav_stats.ROL_peaks);

idx_asian = (HFB_behav_stats.test_race ==1);
idx_black = (HFB_behav_stats.test_race ==2);
idx_white = (HFB_behav_stats.test_race ==3);

data_asian = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_asian,HFB_range(1):HFB_range(2));
data_black = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_black,HFB_range(1):HFB_range(2));
data_white = HFB_behav_stats.raw_HFB(idx_clean&idx_ROL_nan&idx_white,HFB_range(1):HFB_range(2));

data_asian_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_asian,:);
[data_asian_ROL_peaks_s,asian_sort] = sort(data_asian_ROL_peaks);
data_asian_s = data_asian(asian_sort,:);
asian_peaks_onsets = HFB_behav_stats.ROL_onsets(asian_sort,:);

data_black_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_black,:);
[data_black_ROL_peaks_s,black_sort] = sort(data_black_ROL_peaks);
data_black_s = data_black(black_sort,:);
black_peaks_onsets = HFB_behav_stats.ROL_onsets(black_sort,:);

data_white_ROL_peaks = HFB_behav_stats.ROL_peaks(idx_clean&idx_ROL_nan&idx_white,:);
[data_white_ROL_peaks_s,white_sort] = sort(data_white_ROL_peaks);
data_white_s = data_white(white_sort,:);
white_peaks_onsets = HFB_behav_stats.ROL_onsets(white_sort,:);

% plot results

figure('Position', [1500 500 450 600]),clf

subplot(3,1,1)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_asian,1):-1:1,data_asian_s,100,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_asian_ROL_peaks_s,size(data_asian,1):-1:1,'k--','LineWidth',3)
plot(asian_peaks_onsets,size(data_asian,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Asian'])

subplot(3,1,2)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_black,1):-1:1,data_black_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_black_ROL_peaks_s,size(data_black,1):-1:1,'k--','LineWidth',3)
plot(black_peaks_onsets,size(data_black,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, Black'])

subplot(3,1,3)
contourf(data_all.time(HFB_range(1):HFB_range(2)),size(data_white,1):-1:1,data_white_s,40,'linecolor','none')
set(gca,'clim',[-1 1.5])
hold on
plot(data_white_ROL_peaks_s,size(data_white,1):-1:1,'k--','LineWidth',3)
plot(white_peaks_onsets,size(data_white,1):-1:1,'ro','LineWidth',3)
title([ 'Time-Trial-HFB plot of ' num2str(anat_name) ' sorted by ROL peaks, White'])

data_stats = [data_asian_ROL_peaks;data_black_ROL_peaks;data_white_ROL_peaks];
data_groups = [ones(length(data_asian_ROL_peaks),1);ones(length(data_black_ROL_peaks),1)*2;ones(length(data_white_ROL_peaks),1)*3];
[pval,anovatable,stats] = anova1(data_stats,data_groups);
multcompare(stats,'ctype','lsd')
%%
switch anat_name
    case 'INSULA'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,2,5,6],:);
    T3.anat{3} = T3.anat{3}(3);
        else
        end
    case 'IFS IFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([5,10],:);
            T3.anat{1} = T3.anat{1}(5);
        else
        end
    case 'SFS SFG'
        if strcmp(site_pick,'One site  per subj')
            channame = channame([1,4,6],:);
            T3.anat{2} = T3.anat{2}(3);
            T3.anat{3} = T3.anat{3}(1);
        else
        end
end
disp(channame);
%%

