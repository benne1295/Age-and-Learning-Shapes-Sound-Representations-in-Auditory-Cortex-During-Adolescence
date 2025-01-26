clearvars -except EducageTable_cell
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , ['Select behavior: EducageTable_cell ']);
addpath(path)
load (file)
cd(path)

%%
clearvars -except EducageTable_cell
binsize = 100; % binsize for dprime trajectory
dprime_threshold = 1; % trials to threshold
numtrials = 1000; % naive and expert trial number
max_ISI = 60; % maximal inter trial interval
smoothing_window = 20 ; % smoothing of the trajectory

% colors
color_age = {[0.5 0.5 0.5],[0 0 0]}
color_eh = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_eh_alpha = {[0.3010 0.7450 0.9330 0.4 ], [0 0.4470 0.7410 0.4 ]} ;
color_lick = [0.4660 0.6740 0.1880] ;
color_no_lick = [0.6350 0.0780 0.1840] ;
color_lick_back = [0.4660 0.6740 0.1880 0.1] ;
color_no_lick_back = [0.6350 0.0780 0.1840 0.1] ;
color_mean_marker = {[0.5 0.5 0.5 0.8]} ;
color_plot_marker = {[0.5 0.5 0.5 0.2]} ;


% line styles 
L = {['-'], ['-'],['-'], ['-']} ;
L2 = {['--'], ['-'],['--'], ['-']};
M = {['o'],['o'],['o'],['o']} ;


group_string = {'adolescent','adult'}
directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';

%% calculate the number of trials
% The lion's den of dprime 
clc
tic

    Mean_adolescent_dPrime_easy = nan (20,numel(EducageTable_cell)) ;
    Mean_adult_dPrime_easy= nan (20,numel(EducageTable_cell)) ;
    Mean_adolescent_dPrime_hard= nan (20,numel(EducageTable_cell)) ;
    Mean_adult_dPrime_hard= nan (20,numel(EducageTable_cell)) ;

for i   = 1:numel(EducageTable_cell)

    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);


    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);
       % [mouse_table] = max_trial_dist (mouse_table, max_ISI);
        mice_name{n,i} = unique(string(mouse_table.mouse_name));

        easy_length(n,i) = length(find(mouse_table.level == 3)) ;

        hard_easy_length(n,i) = length(find (mouse_table.level == 4 & (mouse_table.stimID ==1 | mouse_table.stimID == 4)));
        hard_hard_length(n,i) = length(find (mouse_table.level == 4 & (mouse_table.stimID ==2  | mouse_table.stimID == 3)));

        mean_easy_i(i) = round(nanmean(easy_length(:,i))) ;
        stde_easy_i(i) = round(nanstd(easy_length(:,i)) / sqrt(size(easy_length(:,i),1))) ;

        mean_hard_i(i) = round(nanmean(hard_hard_length(:,i))) ;
        stde_hard_i(i) = round(nanstd(hard_hard_length(:,i)) / sqrt(size(hard_hard_length(:,i),1))) ;

    end
 

    easy_length(easy_length == 0) = nan ;
    hard_easy_length(hard_easy_length == 0) = nan ;
    hard_hard_length(hard_hard_length == 0) = nan ;
end 
toc
 %% calculate dprime under different criteria
clc
tic
for i   = 1:numel(EducageTable_cell)

    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)

        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        %exclude trials above max ISI
        %[mouse_table] = max_trial_dist (mouse_table, max_ISI);

         cage_ID (n,i) = unique(mouse_table.table_num) ;


        % trial outcomes throughout the task discrimination (1 oct. and 0.25 oct. together)
        index_all = find (mouse_table.level == 3 | mouse_table.level == 4);
        [hit_all{n,i}, fa_all{n,i},~, cr_all{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        % trial outcomes in the easy level:  1 oct. discrimination
        index_easy_level =find (mouse_table.level == 3);
        [hit_easy_level{n,i}, fa_easy_level{n,i}, ~, cr_easy_level{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy_level, mouse_table);
        [~, ~, dPrime_easy_level{n,i}, cbias_easy_level{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);
        
        % trial outcomes in the hard level + 1 oct. + 0.25 oct. discrimination
        index_easy = find (mouse_table.level == 4 & (mouse_table.stimID ==1 | mouse_table.stimID == 4));
        index_hard = find (mouse_table.level == 4 & (mouse_table.stimID ==2  | mouse_table.stimID == 3));

        % 1 oc.t discrimination
        [hit_easy{n,i}, fa_easy{n,i}, ~, cr_easy{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy, mouse_table);
        [~, ~, dPrime_easy{n,i}, cbias_easy{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);
        % 0.25 oct. discrimination
        [hit_hard{n,i}, fa_hard{n,i}, ~, cr_hard{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_hard, mouse_table);
        [~, ~, dPrime_hard{n,i}, cbias_hard{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        % trial to threshold - procedural learning
        trial_criterion_easy_level = find((dPrime_easy_level{n,i}(2:end-1)>= dprime_threshold))+1;
        if ~isempty(trial_criterion_easy_level);
            percentage_trial_criterion_easy_level(n,i) = length(trial_criterion_easy_level) / length(dPrime_easy_level{n,i});

            trial_criterion_easy_all_level(n,i) = trial_criterion_easy_level(1) ;%/ length(dPrime_easy_level{n,i});

        elseif isempty (trial_criterion_easy_level)
            trial_criterion_easy_all_level(n,i) = NaN;
        end

        dprime_easy_all{n,i} = [dPrime_easy_level{n,i} dPrime_easy{n,i}]  ;

        Naive_dPrime_easy(n,i) = dPrime_easy_level{n,i}(2) ;
        Naive_dPrime_hard(n,i) = dPrime_hard{n,i}(2) ;

        % unbiased dprime
        min_easy_length = round(min(easy_length(:,i))) ; 
        Min_dPrime_easy(n,i) = dPrime_easy_level{n,i}(round(round(min_easy_length)/binsize)-1) ;

        min_hard_length = round(min(hard_hard_length(:,i))) ; 
        Min_dPrime_hard(n,i) = dPrime_hard{n,i}(round(round(min_hard_length)/binsize)-1);
     

        if     length(dprime_easy_all{n,i}) >= round(mean_easy_i(1)/binsize)-1
            Mean_adolescent_dPrime_easy(n,i) = dprime_easy_all{n,i}(round(mean_easy_i(1)/binsize)-1);
        elseif length(dprime_easy_all{n,i}) < round(mean_easy_i(1)/binsize)-1
               Mean_adolescent_dPrime_easy(n,i) = NaN ; 
        end
        if     length(dprime_easy_all{n,i}) >= round(mean_easy_i(2)/binsize)-1
            Mean_adult_dPrime_easy(n,i) = dprime_easy_all{n,i}(round(mean_easy_i(2)/binsize)-1) ; 
        elseif    length(dprime_easy_all{n,i}) < round(mean_easy_i(2)/binsize)-1
   Mean_adult_dPrime_easy(n,i) = NaN ; 
        end
        if length(dPrime_hard{n,i}) >= round(mean_hard_i(1)/binsize)-1
            Mean_adolescent_dPrime_hard(n,i) = dPrime_hard{n,i}(round(mean_hard_i(1)/binsize)-1);
        elseif length(dPrime_hard{n,i}) <  round(mean_hard_i(1)/binsize)-1
             Mean_adolescent_dPrime_hard(n,i) = NaN ; 
        end
        if length(dPrime_hard{n,i}) >= round(mean_hard_i(2)/binsize)-1
            Mean_adult_dPrime_hard(n,i) = dPrime_hard{n,i}(round(mean_hard_i(2)/binsize)-1) ;
        elseif length(dPrime_hard{n,i}) < round(mean_hard_i(2)/binsize)-1
             Mean_adult_dPrime_hard(n,i) = NaN ; 
        end

        Expert_dPrime_easy_level(n,i) = dPrime_easy_level{n,i}(end);
        Expert_dPrime_easy(n,i) = dPrime_easy{n,i}(end);
        Expert_dPrime_hard(n,i) =  dPrime_hard{n,i}(end);

        if Naive_dPrime_easy(n,i) > 0
            Delta_dprime_easy(n,i) = Expert_dPrime_easy(n,i) - Naive_dPrime_easy(n,i) ;
        elseif Naive_dPrime_easy(n,i) <= 0
            Delta_dprime_easy(n,i) = Expert_dPrime_easy(n,i) + Naive_dPrime_easy(n,i) ;
        end

        if Naive_dPrime_hard(n,i) > 0
            Delta_dprime_hard(n,i) = Expert_dPrime_hard(n,i) - Naive_dPrime_hard(n,i) ;
        elseif Naive_dPrime_hard(n,i) <= 0
            Delta_dprime_hard(n,i) = Expert_dPrime_hard(n,i) + Naive_dPrime_hard(n,i) ;
        end


        if Naive_dPrime_easy(n,i) > 0
            Delta_min_dprime_easy(n,i) = Min_dPrime_easy(n,i) - Naive_dPrime_easy(n,i) ;
        elseif Naive_dPrime_easy(n,i) <= 0
            Delta_min_dprime_easy(n,i) = Min_dPrime_easy(n,i) + Naive_dPrime_easy(n,i) ;
        end

        if Naive_dPrime_hard(n,i) > 0
            Delta_min_dprime_hard(n,i) = Min_dPrime_hard(n,i) - Naive_dPrime_hard(n,i) ;
        elseif Naive_dPrime_hard(n,i) <= 0
            Delta_min_dprime_hard(n,i) = Min_dPrime_hard(n,i) + Naive_dPrime_hard(n,i) ;
        end

        Delta_switch(n,i) = dPrime_easy{n,i}(1) - dPrime_easy_level{n,i}(end) ; 
        Delta_switch_cbias(n,i) = cbias_easy{n,i}(1) - cbias_easy_level{n,i}(end) ; 

    end

end
toc 
 %%  Figure 1 D: learning curve of an example animal 
clc
close all
% =i 1 n = 12
%i = 2 n = 5
for i   = 1:numel(EducageTable_cell)
 %fig1d = figure
  EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);

    for n= 1:length(mice)
%subplot(4,5,n)
figure
    nan_easy_level = [] ;
    d_easy = [dPrime_easy_level{n,i} dPrime_easy{n,i}]  ;
    nan_easy_level = nan(1,size(dPrime_easy_level{n,i},2)) ;
    d_hard = [nan_easy_level dPrime_hard{n,i}] ;

    plot(d_easy,'Color',color_eh{1},'linestyle',L{i},'linewidth',1)
    hold on
    plot(d_hard,'Color',color_eh{2},'linestyle',L{i},'linewidth',1)
    hold on
    ylim([-1 5])
    xlim([0 150])
    xline(size(nan_easy_level,2),'--')
    yline(1,'--')
    xticks([0:50:150])
    xticklabels ([0:5000:15000])
    yticks([-1:2:5])
    yticklabels ([-1:2:5])
    xlabel('trials')
    ylabel("d'")
    box off;
    hold on
    ax = gca;
    %ax.XAxis.FontSize = 20;
    %ax.YAxis.FontSize = 20;
    movegui('east');
    title (mice_name{n,i})
    end


end 

    filename = 'Fig1C_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(7);
saveas(recent_figure, filepath, 'svg');
   filename = 'Fig1C_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(25);
saveas(recent_figure, filepath, 'svg');

%% trials to threshold d' = 1
close all
clc
f1e = figure
trial_criterion_easy_all_level(trial_criterion_easy_all_level == 0) = nan ;
violinplot([ trial_criterion_easy_all_level(:,1),trial_criterion_easy_all_level(:,2)],{'adol.','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} ) ;
xticks ([1 2])
xlim([0 3]);
%ylim([0 1.2]);
%yticks([0:2000:6000])
%yticklabels([0:2000:6000])
xticklabels({'adolescent','adult'})
ylabel(" trials to criterion")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


[p,h,stats] = ranksum(trial_criterion_easy_all_level(:,1), trial_criterion_easy_all_level(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95] ;
line_height = [0.8] ;
text_mean = [1.5] ;
text_height = [1]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 
line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

   end

    filename = 'Fig1D'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% overview of level transition, and mean dprime at level transtiion 
clc
close all
for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);
    for n=1:length(mice)
        dprime_level_end_easy(n,i) = dPrime_easy_level{n,i}(end) ;
       n_trials(n,i) =  length(dPrime_easy_level{n,i}) ;
    end

    mean_dprime_end_easy(i) = nanmean(dprime_level_end_easy(:,i)) ;
    stde_dprime_end_easy(i) = nanstd(dprime_level_end_easy(:,i))/ sqrt(numel(dprime_level_end_easy(:,i))); 

    for n=1:length(mice)
        figure(1)
        nan_easy_level = [] ;
        plot(length(dPrime_easy_level{n,i}),dPrime_easy_level{n,i}(end),'Color',color_age{i},'linestyle','none','Marker','o'...
            ,'MarkerSize',10,'MarkerEdgeColor',color_age{i}, 'MarkerFaceColor',color_age{i})
        hold on
        ylim([-1 5])
        xlim([0 200])
        %xline(length(dPrime_easy_level{n,i}),'--','linewidth',1,'Color',color_age{i})
        xline(round(mean_easy_i(i)/100),'--','linewidth',3,'Color',color_age{i})
        yline( mean_dprime_end_easy(i) ,'--','linewidth',3,'Color',color_age{i})
        yline(1,'--','linewidth',1)
        xticks([0:50:200])
        xticklabels ([0:5000:20000])
        yticks([-1:2:5])
        yticklabels ([-1:2:5])
       xlabel('#trials at P30 ')
        ylabel("d' at P30")
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
    end
end

filename = 'Fig1E'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

[h,p,ci,stats] = ttest2(n_trials(:,1), n_trials(:,2)) ;

[p,h,stats] = ranksum(dprime_level_end_easy(:,1), dprime_level_end_easy(:,2)) ;

[R,P] = corrcoef(n_trials(:,1),dprime_level_end_easy(:,1))
[R,P] = corrcoef(n_trials(:,2),dprime_level_end_easy(:,2))

%% minimal trials violinplot easy level
close all
clc
figure
violinplot([Min_dPrime_easy(:,1),....
  Min_dPrime_easy(:,2)],...
   {group_string{1},group_string{2}},"ViolinColor",{[color_eh{1};color_eh{1}]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-1 5]);
yticks([-1 1 3  5])
        yline(1,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(Min_dPrime_easy(:,1), Min_dPrime_easy(:,2)) ;

line_left = [1.05 ] ;
line_right =  [1.95 ] ;
line_height = [4.3 ] ;
text_mean = [1.5] ;
text_height = [4.5]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
   end

filename = 'sUPPFig11c'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');


%% mean adolescent trials violinplot easy level
close all
clc
figure
violinplot([Mean_adolescent_dPrime_easy(:,1),....
    Mean_adolescent_dPrime_easy(:,2)],...
    {group_string{1},group_string{2}},"ViolinColor",{[color_eh{1};color_eh{1}]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-2 6]);
yticks([-2 -1 0 1 2 3 4 5])
yline(1,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(Mean_adolescent_dPrime_easy(:,1), Mean_adolescent_dPrime_easy(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [4.3] ;
text_mean = [1.5 ] ;
text_height = [4.5 ]  ;

s1 = [1] ;
for s = 1:length(s1)
    txt = num2str(p(s1(s))) ;

    line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
    hold on
    text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
end

   filename = 'SuppFig11d'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% mean adult trials violinplot easy level
close all
clc
figure
violinplot([Mean_adult_dPrime_easy(:,1),....
    Mean_adult_dPrime_easy(:,2)],...
    {group_string{1},group_string{2}},"ViolinColor",{[color_eh{1};color_eh{1}]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-2 6]);
yticks([-2 -1 0 1 2 3 4 5])
yline(1,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(Mean_adult_dPrime_easy(:,1), Mean_adult_dPrime_easy(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [4.3] ;
text_mean = [1.5 ] ;
text_height = [4.5 ]  ;

s1 = [1] ;
for s = 1:length(s1)
    txt = num2str(p(s1(s))) ;

    line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
    hold on
    text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
end

   filename = 'SuppFig11e'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');



  %% violionplot delta dprime between dprime of easy task after level switch  
   close all
clc
figure
violinplot([Delta_switch(:,1),....
  Delta_switch(:,2)],...
   {group_string{1},group_string{2}},"ViolinColor",{[color_eh{1};color_eh{1}]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-3 5]);
yticks([-3 -2 -1 0 1 2 3 ])
      %  yline(1,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("easy \Delta d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(Delta_switch(:,1), Delta_switch(:,2)) ;

[p1,h,stats] = signrank(Delta_switch(:,1)) ;

[p2,h,stats] = signrank(Delta_switch(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [3.3] ;
text_mean = [1.5 ] ;
text_height = [3.5]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 
                    txt1 = num2str(p1(s1(s))) ; 
          txt2 = num2str(p2(s1(s))) ; 


line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

text(mean([line_left(s) line_left(s)]), text_height(s) ,txt1,'FontSize',5)

text(mean([line_right(s) line_right(s)]), text_height(s) ,txt2,'FontSize',5)

   end

   filename = 'Fig1f'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% overview of number of trials in the hard task
clc
close all
for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);
    for n=1:length(mice)
        dprime_easy_end(n,i) = dPrime_easy{n,i}(end) ;
        dprime_hard_end(n,i) = dPrime_hard{n,i}(end) ;
       n_trials(n,i) =  length(dPrime_hard{n,i}) ;
    end 

    mean_dprime_end_hard(i) = nanmean(dprime_hard_end(:,i)) ;
    stde_dprime_end_hard(i) = nanstd(dprime_hard_end(:,i))/ sqrt(numel(dprime_hard_end(:,i))); 

    for n=1:length(mice)
        figure(1)
        nan_easy_level = [] ;
        plot(length(dPrime_hard{n,i}),dPrime_hard{n,i}(end),'Color',color_age{i},'linestyle','none','Marker','o'...
            ,'MarkerSize',10,'MarkerEdgeColor',color_age{i}, 'MarkerFaceColor',color_age{i})
        hold on
        ylim([-1 5])
        xlim([0 100])
        %xline(length(dPrime_easy_level{n,i}),'--','linewidth',1,'Color',color_age{i})
        xline(round(mean_hard_i(i)/100),'--','linewidth',3,'Color',color_age{i})
        yline( mean_dprime_end_hard(i) ,'--','linewidth',3,'Color',color_age{i})
        yline(1,'--','linewidth',1)
        xticks([0:25:100])
        xticklabels ([0:2500:10000])
        yticks([-1:2:5])
        yticklabels ([-1:2:5])
       xlabel('#trials of the hard task at P38 ')
        ylabel("d' at P38")
        box off;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
    end
end

[h,p,ci,stats] = ttest2(n_trials(:,1), n_trials(:,2)) ;

[p,h,stats] = ranksum(dprime_easy_end(:,1), dprime_easy_end(:,2)) ;

[p,h,stats] = ranksum(dprime_hard_end(:,1), dprime_hard_end(:,2)) ;

[R,P] = corrcoef(n_trials(:,1),dprime_hard_end(:,1))
[R,P] = corrcoef(n_trials(:,2),dprime_hard_end(:,2))
   

filename = 'Fig1g'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');


   
   %% delta dprime easy violinplot

Delta_dprime_easy(Delta_dprime_easy == 0 ) = nan ; 

close all
clc
figure
violinplot([Delta_dprime_easy(:,1),....
   Delta_dprime_easy(:,2)],...
   {group_string{1},group_string{2}},"ViolinColor",{[color_eh{1};color_eh{1}]}) ;
xticks ([1 2])
xlim([0 3]);

ylim([-2 6]);
yticks([-2 -1 0 1 2 3 4 5])
        yline(0,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("\Delta d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
[p,h,stats] = ranksum(Delta_dprime_easy(:,1), Delta_dprime_easy(:,2)) ;

[p1,h,stats] = signrank(Delta_dprime_easy(:,1)) ;

[p2,h,stats] = signrank(Delta_dprime_easy(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [3.3] ;
text_mean = [1.5 ] ;
text_height = [3.5]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 
                    txt1 = num2str(p1(s1(s))) ; 
          txt2 = num2str(p2(s1(s))) ; 


line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

text(mean([line_left(s) line_left(s)]), text_height(s) ,txt1,'FontSize',5)

text(mean([line_right(s) line_right(s)]), text_height(s) ,txt2,'FontSize',5)

   end

   filename = 'Fig1g'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');


    %% delta dprime hard violinplot

   Delta_dprime_hard(Delta_dprime_hard == 0 ) = nan ; 

close all
clc
figure
violinplot([Delta_dprime_hard(:,1),....
   Delta_dprime_hard(:,2)],...
   {group_string{1},group_string{2}},"ViolinColor",{[color_eh{2};color_eh{2}]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-2 6]);
yticks([-2 -1 0 1 2 3 4 5])
        yline(0,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel("\Delta d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
[p,h,stats] = ranksum(Delta_dprime_hard(:,1), Delta_dprime_hard(:,2)) ;

[p1,h,stats] = signrank(Delta_dprime_hard(:,1)) ;

[p2,h,stats] = signrank(Delta_dprime_hard(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [3.3] ;
text_mean = [1.5 ] ;
text_height = [3.5]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 
                    txt1 = num2str(p1(s1(s))) ; 
          txt2 = num2str(p2(s1(s))) ; 


line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

text(mean([line_left(s) line_left(s)]), text_height(s) ,txt1,'FontSize',5)

text(mean([line_right(s) line_right(s)]), text_height(s) ,txt2,'FontSize',5)

   end

   filename = 'Fig1h'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% Table 1 - LME effect of age,  group housing and sex on dprime
clc
for i   = 1:numel(EducageTable_cell) % perform per group

    dprime_average(:,i) = mean ([Expert_dPrime_easy(:,i)'; Expert_dPrime_hard(:,i)']) ;

    idx_f = find(EducageTable_cell{1,i}.sex == 'female') ;
    idx_m = find(EducageTable_cell{1,i}.sex == 'male') ;

    female_mice = unique(EducageTable_cell{1,i}.mouse_num(idx_f,:)) ;
    male_mice =unique(EducageTable_cell{1,i}.mouse_num(idx_m,:)) ;

    male_ID = nan(height(female_mice),1) ;
    female_ID = nan(height(male_mice),1) ;
    male_ID(1:end,:) = 2 ;
    female_ID(1:end,:) = 1 ;

    age_sex_ID(:,i) = [male_ID' female_ID']
end

%LME test if the effect is relevent per group
adult_ID = nan(height(Expert_dPrime_easy),1) ;
adolescent_ID = nan(height(Expert_dPrime_easy),1) ;
adult_ID(1:end,:) = 2 ;
adolescent_ID(1:end,:) = 1 ;

LME_edu = table ;
LME_edu.group = [adolescent_ID' adult_ID']' ;
LME_edu.dprime = [dprime_average(:,1)' dprime_average(:,2)']' ;
LME_edu.sex = [age_sex_ID(:,1)' age_sex_ID(:,2)']' ;
LME_edu.cage = [cage_ID(:,1)' (cage_ID(:,2)' + max(cage_ID(:,1)))]' ;

% Convert categorical variables
LME_edu.group = categorical(LME_edu.group);
LME_edu.sex = categorical(LME_edu.sex);
LME_edu.cage = categorical(LME_edu.cage);


model_edu = fitlme(LME_edu, 'dprime ~ group + sex + (1|cage)') ;
disp(model_edu);

% Extracting necessary information from the LinearMixedModel object
fixedEffects = model_edu.Coefficients;
%post_hoc comparison
fixedEffects.pValue = fixedEffects.pValue*2
fixedEffects.pValue(fixedEffects.pValue >1) = 0.9999 ; 
fixedEffects.Estimate = round(fixedEffects.Estimate,4)
fixedEffects.SE = round(fixedEffects.SE,4)
fixedEffects.tStat = round(fixedEffects.tStat,4)
fixedEffects.pValue = round(fixedEffects.pValue,4)
fixedEffects.Lower = round(fixedEffects.Lower,4)
fixedEffects.Upper = round(fixedEffects.Upper,4)
randomEffects = model_edu.randomEffects;

% Create a table for fixed effects
fixedEffectsTable = table(fixedEffects.Estimate, fixedEffects.SE, fixedEffects.tStat, fixedEffects.DF, fixedEffects.pValue,  fixedEffects.Lower,...
    fixedEffects.Upper,...
    'VariableNames', {'Estimate', 'StandardError', 'tStatistic', 'DF', 'pValue', 'CI lower', 'CI upper'}, ...
    'RowNames', fixedEffects.Name);

writetable(fixedEffectsTable,'LME_edu.xlsx','Sheet',1)

%% Save for later 
for i = 1:numel(EducageTable_cell)
    % sort out mice with d' easy <1 to ensure the same criterion as in the
    % recording
    dprime_edu.dprime_easy(1:size(Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i),1),i)=  Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i);
    dprime_edu.dprime_hard(1:size(Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i),1),i) =  Expert_dPrime_hard(Expert_dPrime_easy(:,i)>1,i);
end
save('Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\data\dprime_edu.mat', 'dprime_edu')


%% Educage Control table
%% supplementary Figure 1.1 :  control experiment - hard task only
clc
close all
% choose EducageTable_control_cell.mat
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\Educage\EducageTable_cell.mat'...
    , ['Select table: Educage_control_cell ']);
addpath(path)
load (file)
%% calculate dprime in the control task
tic
for i   = 1:numel(EducageTable_control_cell)
    EducageTable = cell2table(EducageTable_control_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        index_hard = find (mouse_table.level == 3 & (mouse_table.stimID ==1  | mouse_table.stimID == 2));

        trial_length(n,i) = length(index_hard) ;

        [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_hard, mouse_table);
        [~, ~, dPrime_hard_control, ~] = dprime(binsize, go_licks, ngo_licks, mice);

        dPrime_control_all{n,i} = dPrime_hard_control;

        Naive_dPrime_control(n,i) = dPrime_hard_control(1) ;

        Expert_dPrime_control(n,i) = dPrime_hard_control(end);

        Delta_Prime_control(n,i) = dPrime_hard_control(end) - dPrime_hard_control(1) ;

        data_plot_control{1,i}(n,:) =  [Naive_dPrime_control(n,i) Expert_dPrime_control(n,i)];

    end
end
toc

%% supplementary Figure 1.1 a: d' control trajectory

for i   = 1:numel(EducageTable_control_cell)
    n = max(cellfun(@length,dPrime_control_all));

    figure
    for n=1:max(mice)

        d_hard = [dPrime_control_all{n,i}] ;

        d_hard = smoothdata(d_hard,'gaussian',smoothing_window)

        plot(d_hard,'Color',color_eh{2},'linestyle',L{i},'linewidth',4)
        hold on
        ylim([-1 5])
        xlim([0 101])
        yline(1,'--')
        xticks([0:50:101])
        xticklabels ([0:5000:10000])
        yticks([-1:2:5])
        yticklabels ([-1:2:5])
        xlabel('trials')
        ylabel("d'")
        box off;
        makepretty;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');


    end
end


   filename = 'Supp_Fig11_a'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% supplementary Figure 1.1 B: d' control naive a expert d'
clc
close all

figure
plot([1 2],nanmean(data_plot_control{1,1}),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5) %adolescent
hold on

plot([1 2],data_plot_control{1,1},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adolescent
hold on

plot([3 4],nanmean(data_plot_control{1,2}),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5) %adolescent
hold on

plot([3 4],data_plot_control{1,2},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adult
hold on
xticks ([1 2 3 4])
xlim([0 5]);
yticks ([-1:5])
ylim([-1 5]);
ylabel("d'")
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.',})

box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
h = kstest(data_plot_control{1,1})
if h > 0
    [p_hard(3,1),~,~] = ranksum(data_plot_control{1,1}(:,1),data_plot_control{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [p_hard(4,1),~,~] = ranksum(data_plot_control{1,1}(:,2),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(3:4,1),~]= bonferroni_holm(p_hard(3:4,1),0.05) ;
    [p_hard(1,1),~,~] = signrank(data_plot_control{1,1}(:,1),data_plot_control{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [p_hard(2,1),~,~] = signrank(data_plot_control{1,2}(:,1),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;

end
disp(p_hard)

if  p_hard(1,1) < 0.0005
    txt = '***';
elseif  p_hard(1,1) < 0.005
    txt = '**';
elseif  p_hard(1,1) < 0.05
    txt = '*';
elseif p_hard(1,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(2,1) < 0.0005
    txt = '***';
elseif  p_hard(2,1) < 0.005
    txt = '**';
elseif  p_hard(2,1) < 0.05
    txt = '*';
elseif p_hard(2,1) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(3,1) < 0.0005
    txt = '***';
elseif  p_hard(3,1) < 0.005
    txt = '**';
elseif  p_hard(3,1) < 0.05
    txt = '*';
elseif p_hard(3,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p_hard(4,1) < 0.0005
    txt = '***';
elseif  p_hard(4,1) < 0.005
    txt = '**';
elseif  p_hard(4,1) < 0.05
    txt = '*';
elseif p_hard(4,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

   filename = 'Supp_Fig11_b'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% psychometric curve calculation

[file, path] = uigetfile('*.xlsx',['Select stimulus file (Disc_hard_Educage) ']);
t = readtable(fullfile(path, file)) ;
% determine stimulus IDs for psychometric cruve
level = 4 % ; 
learned_freqs = [1:4];
catch_freqs = [5:11];
stim = [learned_freqs catch_freqs];
freqs =sort(t.Frequency(t.StimulusID(stim)))';
freqs = freqs(2:end-1);
 go_freqs = freqs(1:5) ;
 ngo_freqs = freqs(5:9) ;
length_go_freqs = [1:5] ;
length_ngo_freqs = [5:9] ;
max_ISI = 60 ; 
binsize = 100 


clc
close all
%% calculate lick rates of last 2k
if ~any(strcmp(EducageTable_cell,'1'))
    for i   = 1:numel(EducageTable_cell)
        EducageTable = cell2table(EducageTable_cell(i));
        EducageTable = EducageTable.Var1{1,1};
        data = EducageTable( EducageTable.level == level, :);

        [data] = max_trial_dist (data, max_ISI);

        mice = unique(EducageTable.mouse_num);

        % extract the lick rates per mouse
        [l_all, l_learned, stderr_l_all, mean_l_all, stderr_l_learend, mean_l_learned] = lick_ratio_all (catch_freqs, learned_freqs,mice, stim,data);
        % extract the mean lick rate per group
        [lick_all{i,1}, stderr_lick_all{i,1}, mean_lick_all{i,1}, lick_learned{i,1}, stderr_lick_learned{i,1}, mean_lick_learned{i,1}] = sort_freqs (level, l_all, stderr_l_all, mean_l_all);
        % fit a sigmoid function and normalize the curve
        % go into the function and uncomment the mat2gray i you wanto to
        % normalize between 0 and 2 as presented in figure 2C
        [fitted_curve_all{i,1},mean_fitted_curve{i,1}, stderr_fitted_curve{i,1}] = fit_norm_sigmoid (freqs,lick_all{i,1});

    end
end


%% Fig.2b : plot raw psychometric curves   (also supplementary Fiure 1.2 B)
clc
close all
for i   = 1:numel(EducageTable_cell)

    figure(3);
    plot(go_freqs,nanmean(lick_all{i,1}(:,1:5)),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
    hold on
    plot(ngo_freqs,nanmean(lick_all{i,1}(:,5:9)),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
    hold on

    patch([go_freqs flip(go_freqs)] , [(nanmean(lick_all{i,1}(:,1:5))) + (stderr_lick_all{i,1}(:,1:5))...
        flip((nanmean(lick_all{i,1}(:,1:5))) - (stderr_lick_all{i,1}(:,1:5)))],color_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    patch([ngo_freqs flip(ngo_freqs)] , [(nanmean(lick_all{i,1}(:,5:9))) + (stderr_lick_all{i,1}(:,5:9))...
        flip((nanmean(lick_all{i,1}(:,5:9))) - (stderr_lick_all{i,1}(:,5:9)))],color_no_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    ylim([0 1]);
    yticks([0:0.2:1]);
    xline(10,'--','Color','k','linewidth',1);
    yline(0.5,'--','Color','k','linewidth',1);

    xlim([min(t.Frequency) max(t.Frequency)])
    xlabel('freq (kHz)');
    ylabel('lick rate');
    set(gca, 'XDir','reverse','XScale','log');
    box off;
    makepretty;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

end

%statistics
h = kstest(lick_all{1,1}) ;
[p_freq(1,1),~,~] = ranksum(reshape(lick_all{1,1}(:,[1 3]),[],1),reshape(lick_all{2,1}(:,[1 3]),[],1),'alpha',0.05,'tail','both')
[p_freq(2,1),~,~] = ranksum(reshape(lick_all{1,1}(:,[7 9]),[],1),reshape(lick_all{2,1}(:,[7 9]),[],1),'alpha',0.05,'tail','both')


[~,p_freq,~]= bonferroni_holm(p_freq,0.05) ;
disp(p_freq)


   filename = 'Fig2b'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% Fig. 2c: normalized psychometric curve fitted to a sigmoid curve
clc
close all

for i   = 1:numel(EducageTable_cell)
    for n = 1:size(fitted_curve_all{i,1},1);
        perc = find (fitted_curve_all{i,1}(n,:) > 0.5);
        percile(n,i) = perc(end);
        percentile(n,i) = fitted_curve_all{1,1}(n,percile(n,i));
        freq_percentile(n,i) = freqs(:,percile(n,i));
    end

    figure(1)
    plot(go_freqs,mean_fitted_curve{i,1}(:,1:5),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
    hold on
    plot(ngo_freqs,mean_fitted_curve{i,1}(:,5:9),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
    hold on

    patch([go_freqs flip(go_freqs)] , [mean_fitted_curve{i,1}(:,1:5) + (stderr_fitted_curve{i,1}(:,1:5))...
        flip(mean_fitted_curve{i,1}(:,1:5) - (stderr_fitted_curve{i,1}(:,1:5)))],color_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    patch([ngo_freqs flip(ngo_freqs)] , [mean_fitted_curve{i,1}(:,5:9) + (stderr_fitted_curve{i,1}(:,5:9))...
        flip(mean_fitted_curve{i,1}(:,5:9) - (stderr_fitted_curve{i,1}(:,5:9)))],color_no_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    ylim([0 1]);
    yticks([0:0.2:1]);
    xline(10,'--','Color','k');
    yline(0.5,'--','Color','k');
    xlim([min(t.Frequency) max(t.Frequency)])
    set(gca, 'XDir','reverse','XScale','log');
    xlabel('freq  (kHz)');
    ylabel('normalized lick rate');
    box off;
    makepretty;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
end

%statistics
[p_freq_perc(1,1),~,~] = ranksum(freq_percentile(:,1),freq_percentile(:,2),'alpha',0.05,'tail','both')
[p_freq_perc(2,1),~,~] = ranksum(percentile(:,1),percentile(:,2),'alpha',0.05,'tail','both')
[~,p_freq_perc,~]= bonferroni_holm(p_freq_perc,0.05) ;
disp(p_freq_perc)


   filename = 'Fig2c'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% Fig 2d: Calculate maximal c-bias
clc
tic

for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)

        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        numtrials_b = int16(length(mouse_table.level(mouse_table.level == 4))/3);
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        index_all = find ( mouse_table.level == 4);
        [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        Expert_bias_all{1,i}(n,:) =  min(cbias_all{n,i}(size(cbias_all{n,i},2):-1:size(cbias_all{n,i},2) -(numtrials_b/binsize)));
    end
end
toc

figure
violinplot([ Expert_bias_all{1,1},Expert_bias_all{1,2}],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} )
xticks ([1 2])
xlim([0 3]);
ylim([-3 0]);
yticks([-3:1:0])
yticklabels([-3:1:0])
xticklabels({'adolescent','adult'})
ylabel(' max. lick bias (c-bias)')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


[p,~,~] = ranksum(Expert_bias_all{1,1},Expert_bias_all{1,2},'alpha',0.05,'tail','both') ;

disp(p)

if p < 0.0005
    txt = '***';
elseif p < 0.005
    txt = '**';
elseif  p < 0.05
    txt = '*';
elseif p > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(-0.2,-0.2),'color','k')
hold on
text(mean([1.5 1.5]), 0 ,txt,'FontSize',20)

   filename = 'Fig2d'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% lick rate psychometric curve level 3

[file, path] = uigetfile('*.xlsx',['Select stimulus file (Disc_easy_Educage) ']);
t = readtable(fullfile(path, file)) ;

level = 3;
learned_freqs = [1:2];
catch_freqs = [3:9];
stim = [learned_freqs catch_freqs];
freqs =sort(t.Frequency(t.StimulusID(stim)))';
freqs = freqs(2:end-1);
go_freqs = freqs(1:4) ;
ngo_freqs = freqs(4:7) ;
length_go_freqs = [1:4] ;
length_ngo_freqs = [4:7] ;

% split psychometric curves in thirds per level 
for th_r = 1:3
    for i   = 1:numel(EducageTable_cell)
        EducageTable = cell2table(EducageTable_cell(i));
        EducageTable = EducageTable.Var1{1,1};
        data = EducageTable( EducageTable.level == level, :);

       [data] = max_trial_dist (data, max_ISI);

        mice = unique(EducageTable.mouse_num);

        % extract the lick rates per mouse
        for  n=1:length(mice)

            mouse_table = data(data.mouse_num== mice(n),:);

            th = int16(length(mouse_table.level)/3);

            if th_r ==1
                th_range = 1:th*1-1 ;
            elseif th_r == 2
                th_range = th*1:th*2-1
            elseif th_r == 3
                th_range = th*2:th*3-1 ;
            end

            mouse_table = mouse_table(th_range,:);
            [lick_ratio_per_stim] =  stim_lick_ratios (learned_freqs, catch_freqs, mouse_table, stim);
            l_all_3{th_r,i}(mice(n),:) = lick_ratio_per_stim;

            index_all = find ( mouse_table.level == level);
            [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
            [~, ~, dPrime_all_3{n,i,th_r}, cbias_all_3{n,i,th_r}] = dprime(binsize, go_licks, ngo_licks, mice);

            Expert_bias_all_3{th_r,i}(n,:) =  min(cbias_all_3{n,i,th_r});

        end

        l_all_3{th_r,i}(isnan(l_all_3{th_r,i}))= 0;
        [mean_l_all_3{th_r,i}, stderr_l_all_3{th_r,i}] = mean_stderr (l_all_3{th_r,i});

        % extract the mean lick rate per group
        [lick_all_3{th_r,i}, stderr_lick_all_3{th_r,i}, mean_lick_all_3{th_r,i}, lick_learned{th_r,i}...
            , stderr_lick_learned_3{th_r,i}, mean_lick_learned_3{th_r,i}] = sort_freqs (level, l_all_3{th_r,i}, stderr_l_all_3{th_r,i}, mean_l_all_3{th_r,i});
        % fit a sigmoid function and normalize the curve
        [fitted_curve_all_3{th_r,i},mean_fitted_curve_3{th_r,i}, stderr_fitted_curve_3{th_r,i}] = fit_norm_sigmoid (freqs,lick_all_3{th_r,i});
    
    end
end 

  %% Fig. 1 h:  psychometric curve fitted to a sigmoid curve level 3
clc
close all
%close all
for th_r = 1:3

    for i   = 1:numel(EducageTable_cell)

        figure(th_r)
        plot(go_freqs,mean_fitted_curve_3{th_r,i}(:,length_go_freqs),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
        hold on
        plot(ngo_freqs,mean_fitted_curve_3{th_r,i}(:,length_ngo_freqs),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
        hold on

        patch([go_freqs flip(go_freqs)] , [mean_fitted_curve_3{th_r,i}(:,length_go_freqs) + (stderr_fitted_curve_3{th_r,i}(:,length_go_freqs))...
            flip(mean_fitted_curve_3{th_r,i}(:,length_go_freqs) - (stderr_fitted_curve_3{th_r,i}(:,length_go_freqs)))],color_lick,...
            'facealpha' , 0.1, 'EdgeColor','none')

        hold on
        patch([ngo_freqs flip(ngo_freqs)] , [mean_fitted_curve_3{th_r,i}(:,length_ngo_freqs) + (stderr_fitted_curve_3{th_r,i}(:,length_ngo_freqs))...
            flip(mean_fitted_curve_3{th_r,i}(:,length_ngo_freqs) - (stderr_fitted_curve_3{th_r,i}(:,length_ngo_freqs)))],color_no_lick,...
            'facealpha' , 0.1, 'EdgeColor','none')

        hold on
        ylim([0 1]);
        yticks([0:0.2:1]);
        xline(10,'--','Color','k');
        yline(0.5,'--','Color','k');
        xlim([min(go_freqs) max(ngo_freqs)])
        set(gca, 'XDir','reverse','XScale','log');
        xlabel('freq  (kHz)');
        ylabel('fraction to lick threshold');
        box off;
        makepretty;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
    end
end
   filename = 'Fig2A_3_first_third'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
   filename = 'Fig2A_3_last_third'
filepath = fullfile(directory, filename);
recent_figure = figure(3)
saveas(recent_figure, filepath, 'svg');
%% maximal C-bias level 3
close all
clc
for th_r = 1:3
 
figure(th_r+3)
violinplot([ Expert_bias_all_3{th_r,1},Expert_bias_all_3{th_r,2}],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} )
xticks ([1 2])
xlim([0 3]);
ylim([-3 0]);
yticks([-3:1:0])
yticklabels([-3:1:0])
xticklabels({'adolescent','adult'})
ylabel(' max. lick bias (c-bias)')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
h = [] ;
[p,~,stats] = ranksum(Expert_bias_all_3{th_r,1},Expert_bias_all_3{th_r,2},'alpha',0.05,'tail','both') %;

txt = num2str(p);

line( linspace(1.05,1.95), linspace(-0.2,-0.2),'color','k')
hold on
text(mean([1.5 1.5]), 0 ,txt,'FontSize',15)
end 
  filename = 'Fig2B_3_first_third'
filepath = fullfile(directory, filename);
recent_figure = figure(4);
saveas(recent_figure, filepath, 'svg');

  filename = 'Fig2B_3_last_third'
filepath = fullfile(directory, filename);
recent_figure = figure(6);
saveas(recent_figure, filepath, 'svg');
%% lick rate psychometric curve level 4

[file, path] = uigetfile('*.xlsx',['Select stimulus file (Disc_hard_Educage) ']);
t = readtable(fullfile(path, file)) ;
% determine stimulus IDs for psychometric cruve
level = 4 % ; 
learned_freqs = [1:4];
catch_freqs = [5:11];
stim = [learned_freqs catch_freqs];
freqs =sort(t.Frequency(t.StimulusID(stim)))';
freqs = freqs(2:end-1);
 go_freqs = freqs(1:5) ;
 ngo_freqs = freqs(5:9) ;
length_go_freqs = [1:5] ;
length_ngo_freqs = [5:9] ;
max_ISI = 60 ; 
binsize = 100 

for th_r = 1:3
    for i   = 1:numel(EducageTable_cell)
        EducageTable = cell2table(EducageTable_cell(i));
        EducageTable = EducageTable.Var1{1,1};
        data = EducageTable( EducageTable.level == level, :);

        [data] = max_trial_dist (data, max_ISI);

        mice = unique(EducageTable.mouse_num);

        % extract the lick rates per mouse
        for  n=1:length(mice)

            mouse_table = data(data.mouse_num== mice(n),:);

            th = int16(length(mouse_table.level)/3);

            if th_r ==1
                th_range = 1:th*1-1 ;
            elseif th_r == 2
                th_range = th*1:th*2-1
            elseif th_r == 3
                th_range = th*2:th*3-1 ;
            end

            mouse_table = mouse_table(th_range,:);
            [lick_ratio_per_stim] =  stim_lick_ratios (learned_freqs, catch_freqs, mouse_table, stim);
            l_all_4{th_r,i}(mice(n),:) = lick_ratio_per_stim;

            index_all = find ( mouse_table.level == level);
            [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
            [~, ~, dPrime_all_4{n,i,th_r}, cbias_all_4{n,i,th_r}] = dprime(binsize, go_licks, ngo_licks, mice);

            Expert_bias_all_4{th_r,i}(n,:) =  min(cbias_all_4{n,i,th_r});

        end

        l_all_4{th_r,i}(isnan(l_all_4{th_r,i}))= 0;
        [mean_l_all_4{th_r,i}, stderr_l_all_4{th_r,i}] = mean_stderr (l_all_4{th_r,i});

        % extract the mean lick rate per group
        [lick_all_4{th_r,i}, stderr_lick_all_4{th_r,i}, mean_lick_all_4{th_r,i}, lick_learned_4{th_r,i}...
            , stderr_lick_learned_4{th_r,i}, mean_lick_learned_4{th_r,i}] = sort_freqs (level, l_all_4{th_r,i}, stderr_l_all_4{th_r,i}, mean_l_all_4{th_r,i});
        % fit a sigmoid function and normalize the curve
        [fitted_curve_all_4{th_r,i},mean_fitted_curve_4{th_r,i}, stderr_fitted_curve_4{th_r,i}] = fit_norm_sigmoid (freqs,lick_all_4{th_r,i});
    
    end
end 

%% Fig. 1 h:  psychometric curve fitted to a sigmoid curve level 4
clc
close all
for th_r = 1:3

    for i   = 1:numel(EducageTable_cell)

        figure(th_r)
        plot(go_freqs,mean_fitted_curve_4{th_r,i}(:,length_go_freqs),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
        hold on
        plot(ngo_freqs,mean_fitted_curve_4{th_r,i}(:,length_ngo_freqs),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
        hold on

        patch([go_freqs flip(go_freqs)] , [mean_fitted_curve_4{th_r,i}(:,length_go_freqs) + (stderr_fitted_curve_4{th_r,i}(:,length_go_freqs))...
            flip(mean_fitted_curve_4{th_r,i}(:,length_go_freqs) - (stderr_fitted_curve_4{th_r,i}(:,length_go_freqs)))],color_lick,...
            'facealpha' , 0.1, 'EdgeColor','none')

        hold on
        patch([ngo_freqs flip(ngo_freqs)] , [mean_fitted_curve_4{th_r,i}(:,length_ngo_freqs) + (stderr_fitted_curve_4{th_r,i}(:,length_ngo_freqs))...
            flip(mean_fitted_curve_4{th_r,i}(:,length_ngo_freqs) - (stderr_fitted_curve_4{th_r,i}(:,length_ngo_freqs)))],color_no_lick,...
            'facealpha' , 0.1, 'EdgeColor','none')

        hold on
        ylim([0 1]);
        yticks([0:0.2:1]);
        xline(10,'--','Color','k');
        yline(0.5,'--','Color','k');
        xlim([min(go_freqs) max(ngo_freqs)])
        set(gca, 'XDir','reverse','XScale','log');
        xlabel('freq  (kHz)');
        ylabel('fraction to lick threshold');
        box off;
        makepretty;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
    end
end
   filename = 'Fig2A_4_first_third'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
   filename = 'Fig2A_4_last_third'
filepath = fullfile(directory, filename);
recent_figure = figure(3)
saveas(recent_figure, filepath, 'svg');
%% maximal C-bias level 4
clc
for th_r = 1:3
 
figure(th_r+3)
violinplot([ Expert_bias_all_4{th_r,1},Expert_bias_all_4{th_r,2}],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} )
xticks ([1 2])
xlim([0 3]);
ylim([-3 0]);
yticks([-3:1:0])
yticklabels([-3:1:0])
xticklabels({'adolescent','adult'})
ylabel(' max. lick bias (c-bias)')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
h = [] ;
[p,~,stats] = ranksum(Expert_bias_all_4{th_r,1},Expert_bias_all_4{th_r,2},'alpha',0.05,'tail','both') %;

txt = num2str(p);

line( linspace(1.05,1.95), linspace(-0.2,-0.2),'color','k')
hold on
text(mean([1.5 1.5]), 0 ,txt,'FontSize',15)
end 

  filename = 'Fig2B_4_first_third'
filepath = fullfile(directory, filename);
recent_figure = figure(4);
saveas(recent_figure, filepath, 'svg');

  filename = 'Fig2B_4_last_third'
filepath = fullfile(directory, filename);
recent_figure = figure(6);
saveas(recent_figure, filepath, 'svg');
   %% violionplot delta cbias between level 3 and level 4
   close all
clc
figure
violinplot([Delta_switch_cbias(:,1),....
  Delta_switch_cbias(:,2)],...
   {group_string{1},group_string{2}},"ViolinColor",{[0.5 0.5 0.5]; [0.5 0.5 0.5]}) ;
xticks ([1 2 ])
xlim([0 3]);

ylim([-2 4]);
      %  yline(1,'--','linewidth',1)

xticklabels({group_string{1},group_string{2}})
ylabel(" \Delta c-bias")
yline(0,'--')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(Delta_switch_cbias(:,1), Delta_switch_cbias(:,2)) ;


[p1,h,stats] = signrank(Delta_switch_cbias(:,1)) ;

[p2,h,stats] = signrank(Delta_switch_cbias(:,2)) ;

line_left = [1.05] ;
line_right =  [1.95 ] ;
line_height = [3.3] ;
text_mean = [1.5 ] ;
text_height = [3.5]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 
                    txt1 = num2str(p1(s1(s))) ; 
          txt2 = num2str(p2(s1(s))) ; 


line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

text(mean([line_left(s) line_left(s)]), text_height(s) ,txt1,'FontSize',5)

text(mean([line_right(s) line_right(s)]), text_height(s) ,txt2,'FontSize',5)

   end

   filename = 'Fig2C'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% calculate the CV 
for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);
        numtrials_b = int16(length(mouse_table.level(mouse_table.level == 4))/3);
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        index_easy_level = find (mouse_table.level == 3 & (mouse_table.stimID ==1 | mouse_table.stimID == 2));
        index_easy = find (mouse_table.level == 4 & (mouse_table.stimID ==1 | mouse_table.stimID == 4));
        index_hard = find (mouse_table.level == 4 & (mouse_table.stimID ==2  | mouse_table.stimID == 3));

           [~,  ~, ~,cr_3{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy_level, mouse_table);
        [~, ~, dPrime_3{n,i}, cbias_3{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);


           [~,  ~, ~,cr_4_easy{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy, mouse_table);
        [~, ~, dPrime_4_easy{n,i}, cbias_4_easy{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);


           [~,  ~, ~,cr_4_hard{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_hard, mouse_table);
        [~, ~, dPrime_4_hard{n,i}, cbias_4_hard{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        %calculate the CV
        cv_cr_3 (n,i) = abs( std( cr_3{n,i}) / mean( cr_3{n,i}) );
        cv_dprime_3 (n,i) = abs(std( dPrime_3{n,i}) / mean( dPrime_3{n,i})) ;
        cv_cbias_3 (n,i) =  abs(std( cbias_3{n,i}) / mean( cbias_3{n,i})) ;

        cv_cr_4_easy (n,i) = abs( std( cr_4_easy{n,i}) / mean( cr_4_easy{n,i}) );
        cv_dprime_4_easy (n,i) =  abs(std( dPrime_4_easy{n,i}) / mean( dPrime_4_easy{n,i})) ;
        cv_cbias_4_easy (n,i) =  abs(std( cbias_4_easy{n,i}) / mean( cbias_4_easy{n,i})) ;

        cv_cr_4_hard (n,i) = abs( std( cr_4_hard{n,i}) / mean( cr_4_hard{n,i}) );
        cv_dprime_4_hard (n,i) =  abs(std( dPrime_4_hard{n,i}) / mean( dPrime_4_hard{n,i})) ;
        cv_cbias_4_hard (n,i) =  abs(std( cbias_4_hard{n,i}) / mean( cbias_4_hard{n,i})) ;



    end
end

%% CV of dprime during easy level 
figure
violinplot([ cv_dprime_3(:,1),cv_dprime_3(:,2)],{'adol.','adult'},"ViolinColor",{color_eh{1}} )

xticks ([1 2])
xlim([0 3]);
ylim([0 3]);
yticks([0:1:3])
yticklabels([0:1:3])
xticklabels({'adol.','adult'})
ylabel("CV of d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(cv_dprime_3(:,1), cv_dprime_3(:,2)) ;

line_left = [1.05 ] ;
line_right =  [1.95] ;
line_height = [2.3 ] ;
text_mean = [1.5] ;
text_height = [2.5 ]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
   end

      filename = 'Fig2D'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% CV of dprime during easy stimulus after the level switch  

figure
violinplot([ cv_dprime_4_easy(:,1),cv_dprime_4_easy(:,2)],{'adol.','adult'},"ViolinColor",{color_eh{1}} )

xticks ([1 2])
xlim([0 3]);
ylim([0 3]);
yticks([0:1:3])
yticklabels([0:1:3])
xticklabels({'adol.','adult'})
ylabel("CV of d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(cv_dprime_4_easy(:,1), cv_dprime_4_easy(:,2)) ;

line_left = [1.05 ] ;
line_right =  [1.95] ;
line_height = [2.3 ] ;
text_mean = [1.5] ;
text_height = [2.5 ]  ;

s1 = [1] ;
   for s = 1:length(s1)
          txt = num2str(p(s1(s))) ; 

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
   end

         filename = 'Fig2E'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% CV of dprime during hard stimulus after the level switch  


figure
violinplot([ cv_dprime_4_hard(:,1),cv_dprime_4_hard(:,2)],{'adol.','adult'},"ViolinColor",{color_eh{2}} )

xticks ([1 2])
xlim([0 3]);
%ylim([0 3]);
%yticks([0:1:3])
%yticklabels([0:1:3])
xticklabels({'adol.','adult'})
ylabel("CV of d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
[p,h,stats] = ranksum(cv_dprime_4_hard(:,1), cv_dprime_4_hard(:,2)) ;

line_left = [1.05 ] ;
line_right =  [1.95] ;
line_height = [2.3 ] ;
text_mean = [1.5] ;
text_height = [2.5 ]  ;

s1 = [1] ;
for s = 1:length(s1)
    txt = num2str(p(s1(s))) ;

    line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
    hold on
    text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)
end


      filename = 'Fig2F'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');