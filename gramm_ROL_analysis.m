addpath(genpath('/Users/chao/Desktop/function_tools/gramm-master/'))

load cdcol.mat

        
cdcol.asian_red  =             [0.9647 0.1647 0];
cdcol.black_blue =             [0.3020 0.3922 0.5529];
cdcol.white_blue =             [0      0.5800 0.7920];
cdcol.own_race_red =           [0.8200 0      0.1800];
cdcol.other_race_black =       [0      0.1608 0.2353];

plot_params.col = [cdcol.asian_red;
            cdcol.black_blue;
            cdcol.white_blue];
load example_data;
% adjust the sequence of subject
hfb2 = hfb;hfb2(4) = hfb(5);hfb(5) = hfb(4);
race2 = race;race2(4) = race(5);race2(5) = race(4);
sites2 = sites;sites(4) = sites(5);sites(5) = sites(4);
sitescat2 = sitescat;

for i = 1: length(sitescat2)
    if strcmp(sitescat2{i},'S2_X6R')
        sitescat2{i} = 'S3_X6R';
    elseif strcmp(sitescat2{i},'S3_X6R')
        sitescat2{i} = 'S2_X6R';
    else
    end
end

[h1,p1] = ttest2(hfb{i}(strcmp(race{i},'asian')),hfb{i}(strcmp(race{i},'black')));
[h2,p2] = ttest2(hfb{i}(strcmp(race{i},'asian')),hfb{i}(strcmp(race{i},'white')));
[h3,p3] = ttest2(hfb{i}(strcmp(race{i},'black')),hfb{i}(strcmp(race{i},'white')));

figure
clear g
%g=gramm('x',cars.Origin_Region,'y',cars.Horsepower,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
g(1,1) = gramm('x',sitescat2(1:189),'y',hfbcat(1:189),'color',racecat(1:189));
g(2,1) = gramm('x',sitescat2(190:400),'y',hfbcat(190:400),'color',racecat(190:400));

g(1,1).stat_summary('geom',{'bar','black_errorbar'},'setylim',true,'dodge',0.4,'width',0.35);
g(2,1).stat_summary('geom',{'bar','black_errorbar'},'setylim',true,'dodge',0.4,'width',0.35);
g(1,1).no_legend()
g(2,1).no_legend()


g.set_color_options('map',plot_params.col);
g.axe_property('XLabel',[]);
g.axe_property('YLabel',[]);
g.axe_property('Ylim',[0 0.55]);
g.set_text_options('base_size',22,...
    'label_scaling',1,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1);

g.draw();

%% 4:2
figure
clear g
%g=gramm('x',cars.Origin_Region,'y',cars.Horsepower,'color',cars.Cylinders,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
g = gramm('x',sitescat2(1:263),'y',hfbcat(1:263),'color',racecat(1:263));
g = gramm('x',sitescat2(264:400),'y',hfbcat(264:400),'color',racecat(264:400));

g.stat_summary('geom',{'bar','black_errorbar'},'setylim',true,'dodge',0.4,'width',0.35);
g.stat_summary('geom',{'bar','black_errorbar'},'setylim',true,'dodge',0.4,'width',0.35);
g.no_legend()
g.no_legend()


g.set_color_options('map',plot_params.col);
g.axe_property('XLabel',[]);
g.axe_property('YLabel',[]);
g.axe_property('Ylim',[0 0.55]);
g.set_text_options('base_size',22,...
    'label_scaling',1,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1);

g.draw();

%%
%Create variables
figure
clear g
g = gramm('x',racecat,'y',hfbcat,'color',racecat);
g.stat_summary('geom',{'bar','black_errorbar'},'setylim',true,'dodge',1,'width',2);
% g.stat_violin('normalization','width');
g.no_legend()
g.set_color_options('map',plot_params.col);
g.axe_property('XLabel',[]);
g.axe_property('YLabel',[]);
g.axe_property('Ylim',[0 0.55]);

g.set_text_options('base_size',22,...
    'label_scaling',1,...
    'legend_scaling',1,...
    'legend_title_scaling',1,...
    'facet_scaling',1,...
    'title_scaling',1);

g.draw();


% The Stats should be Anova
[p1,tbl1,stats1] = anova1(hfb{i},race{i});
std1 = [std(hfb{i}(strcmp(race{i},'asian'))),std(hfb{i}(strcmp(race{i},'black'))),std(hfb{i}(strcmp(race{i},'white'))) ]
m1 = multcompare(stats1,'ctype','lsd')

[p1,tbl1,stats1] = anova1(hfbcat,racecat);
std1 = [std(hfbcat(strcmp(racecat,'asian'))),std(hfbcat(strcmp(racecat,'black'))),std(hfbcat(strcmp(racecat,'white'))) ]
m1 = multcompare(stats1,'ctype','lsd')
m1 = multcompare(stats1)