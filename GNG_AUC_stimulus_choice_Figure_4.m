clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load GNG_rec_all_cell & Fr_array
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select GNG_rec_all_cell ');
addpath(path)
load (file)

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select Fr_array ');
addpath(path)
load (file)
cd (path)

[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select bin_array ');
addpath(path)
load (file)
cd (path)

%% areas of recording
area_str ={'AUDd','AUDp','AUDv','TEa'};
areas = 1:length(area_str) ;
group_str = {'adolescent','adult'};

% trial types for comparison
% stimulus + choice
% hit vs. fa easy ; hit vs fa hard; cr vs. fa easy; cr vs. fa hard ;
a = [1 7 2 8] ;
b = [2 8 4 10] ;

% AUC parameters
run_window = 50 ; % ms
it_size = 25 ; % ms
n_shuffle = 10 ; % n shuffles for AUC
trial_samples = 10 ;
min_n_trials = 15 ;

%time parameters
window = [-0.2; 0.6]; % time window according to stimulus onset
window_length_ms = dist(window(1),window(2))*1000 ;
n_bins = round( (( (length (1:window_length_ms)) - run_window )/it_size),0) ;

start_base = 1 ;
startStim = 0; % stimulus onset ms
stopStim = 0.1; % stimulus offset  ms
stim_start_ms = dist(window(1), startStim)*1000 ;
stim_stop_ms = dist(window(1), stopStim)*1000 ;
stim_onset_bin = 7 ;


% colors
color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 1]};
color_eh_fade = {[0.3010 0.7450 0.9330 0.3], [0 0.4470 0.7410 0.3], [0 0 1 0.3]};
color_eh_patch = {[0 0.4470 0.7410],[0.3010 0.7450 0.9330], [0 0.4470 0.7410],[0.3010 0.7450 0.9330]};
Colors_area = {[.1 .3 .8], [.5 .4 .9], [0 .5 .6],[0.4940 0.1840 0.5560]};

L = {['--'], ['-'],['--'], ['-']};
M = {['o'],['o'],['o'], ['o']};

directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';


%% only consider experts
GNG_rec_all_cell_exp{1,1} = GNG_rec_all_cell{1,1}; 
GNG_rec_all_cell_exp{1,2} = GNG_rec_all_cell{1,2}; 

%% extract behavior per recording

clc
[behavior] = GNG_neuro_behavior (GNG_rec_all_cell_exp) ;

%% concatenate all cells per recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
tic
[ A_k_ge, B_k_ge, idx_k_ge, T_k_ge, deasy_k_ge, dhard_k_ge, rec_idx, neuron_idx, mouse_idx,area_idx]...
    = GNG_t_per_n_AUC (FR_array, GNG_rec_all_cell_exp, behavior, areas, a, b,  trial_samples, min_n_trials) ;
toc

%% run the AUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[r_AUC, r_AUC_abs, r_AUC_shuf, r_AUC_shuf_abs]...
    = GNG_running_AUC (GNG_rec_all_cell_exp,a,b, A_k_ge, B_k_ge, trial_samples, it_size, run_window, window_length_ms, n_shuffle)
toc
GNG_analysis.r_AUC = r_AUC ;
GNG_analysis.r_AUC_abs = r_AUC_abs ;
GNG_analysis.r_AUC_shuf = r_AUC_shuf ;
GNG_analysis.r_AUC_shuf_abs =  r_AUC_shuf_abs ;


%% mean of stimulus activity  of all trial samples
clc
per_mean = false ;
a = [1 2] ;
stim_onset_bin = abs((window(1)* 1000)/ run_window) ;
[ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell_exp, GNG_analysis.r_AUC_abs, GNG_analysis.r_AUC_shuf_abs,...
    a, trial_samples, run_window, stim_onset_bin, window, it_size, per_mean) ;

GNG_analysis.ind_AUC = ind_AUC ;
GNG_analysis.ind_AUC_shuf  = ind_AUC_shuf ;
GNG_analysis.ind_AUC_all  =ind_AUC_all ;
GNG_analysis.ind_AUC_shuf_all  = ind_AUC_shuf_all ;
GNG_analysis.std_ind_AUC_shuf_all  =std_ind_AUC_shuf_all ;
GNG_analysis.mean_vec =mean_vec ;
GNG_analysis.stde_vec  = stde_vec ;
GNG_analysis.mean_vec_shuf  =mean_vec_shuf ;
GNG_analysis.std_vec_shuf  = std_vec_shuf;

%% mean choice activity of all trial samples
clc
% concatenate e = 3 and e = 4 (fa vs. CR in easy and hard condition
a = [3 4] ; % fa vs. cr mean of easy and hard
per_mean = true; % calculate the mean choice for easy and hard choice together 

[ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell_exp, GNG_analysis.r_AUC_abs, GNG_analysis.r_AUC_shuf_abs,...
    a, trial_samples, run_window, stim_onset_bin, window, it_size, per_mean) ;

for e = a(1)
    for g = 1:numel(GNG_rec_all_cell_exp)
        GNG_analysis.ind_AUC{e,g} = ind_AUC{e,g} ;
        GNG_analysis.ind_AUC_shuf{e,g}  = ind_AUC_shuf{e,g} ;
        GNG_analysis.ind_AUC_all{e,g}  =ind_AUC_all{e,g} ;
        GNG_analysis.ind_AUC_shuf_all{e,g}  = ind_AUC_shuf_all{e,g} ;
        GNG_analysis.std_ind_AUC_shuf_all{e,g}  =std_ind_AUC_shuf_all{e,g} ;
        GNG_analysis.mean_vec{e,g}  =mean_vec{e,g} ;
        GNG_analysis.stde_vec{e,g}   = stde_vec{e,g} ;
        GNG_analysis.mean_vec_shuf{e,g}   =mean_vec_shuf{e,g} ;
        GNG_analysis.std_vec_shuf {e,g}  = std_vec_shuf{e,g};
    end
end

%% Max AUC & AUC Latency Caclulation
clc

a = [1 2 3] ; % stimulus easy ^ hard , choice mean 
[latency_peak_AUC, onset_latency_AUC, width_AUC, max_AUC, rec_idx_th, mouse_idx_th,area_idx_th]...
    = GNG_max_onset_width_AUC (GNG_rec_all_cell_exp, GNG_analysis.ind_AUC_all, GNG_analysis.ind_AUC_shuf_all,...
    GNG_analysis.std_ind_AUC_shuf_all, a, stim_onset_bin,it_size,rec_idx, mouse_idx,area_idx) ;

for e = a
    for g = 1:numel(GNG_rec_all_cell_exp)
        GNG_analysis.max_AUC{e,g} = max_AUC{e,g} ;
        GNG_analysis.width_AUC{e,g}  = width_AUC{e,g} ;
        GNG_analysis.onset_latency_AUC{e,g}  =onset_latency_AUC{e,g} ;
    end
end

 %% Figure 4 E & F  plot a single neuron example 
 % plot hit vs fa easy and hard and choice fa vs cr
 a = [ 1 2 3] ;
 ac = [1 2 3] ; % color code
 close all
 % plot mean (optional include shuffled mean)
 for   g = 2 % expert adult example

     for c = 271 %example neuron
         figure(c)


         for e = 1:length(a)
             subplot(1,3,e)
             plot( GNG_analysis.ind_AUC_all{e,g}(c,:),'Color',color_eh{ac(e)},...
                 'linestyle',L{g},'linewidth',2);
             hold on

             plot( GNG_analysis.ind_AUC_shuf_all{e,g}(c,:),'Color',[0.5 0.5 0.5],...
                 'linestyle',L{g},'linewidth',2)
             hold on
             xline(6,'--','Color','k');
               xticks([0 6 12 18])
             xticklabels([ -200:200:400])
        xlim([1 18])
             ylim([0.45 1])
             yticks([ 0.5 0.75 1])
             xlabel('time (ms)')
             ylabel('Discrimination (abs. AUC)')
             ax = gca;
             ax.XAxis.FontSize = 20;
             ax.YAxis.FontSize = 20;
             movegui('east');
             box off;

         end

         hit_idx = idx_k_ge{2,g}{c}(13,:)  ;
         fa_idx =  idx_k_ge{2,g}{c} (14,:) ;
         miss_idx =  idx_k_ge{2,g}{c} (15,:) ;
         cr_idx  =  idx_k_ge{2,g}{c} (16,:) ;

         hit_idx(hit_idx == 0) =nan ;
         hit_idx =  hit_idx(~isnan(hit_idx)) ;

         fa_idx(fa_idx == 0) =nan ;
         fa_idx =  fa_idx(~isnan(fa_idx)) ;

         miss_idx(miss_idx == 0) =nan ;
         miss_idx =  miss_idx(~isnan(miss_idx)) ;

         cr_idx(cr_idx == 0) =nan ;
         cr_idx =  cr_idx(~isnan(cr_idx)) ;

         color_mean = {[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};

         size_tr = [length(hit_idx),length(miss_idx), length(fa_idx),length(cr_idx)] ;
         max_size = max(size_tr) ;
         trial_types_raster = cat(1,T_k_ge{2,g}{c}(hit_idx,:),T_k_ge{2,g}{c}(miss_idx,:), T_k_ge{2,g}{c}(fa_idx,:),T_k_ge{2,g}{c}(cr_idx,:)) ;

         figure(c+1)
         subplot(2,1,1)
         plot(smooth(mean(T_k_ge{2,g}{c}(hit_idx,:)),20) ,'Color',color_mean{3})
         hold on
         plot(smooth(mean(T_k_ge{2,g}{c}(fa_idx,:)),20) ,'Color',color_mean{2})
         hold on
         plot(smooth(mean(T_k_ge{2,g}{c}(miss_idx,:)),20) ,'Color',color_mean{1})
         hold on
         plot(smooth(mean(T_k_ge{2,g}{c}(cr_idx,:)),20) ,'Color',color_mean{4})
         hold on
         xlabel('time(ms)')
         xline(200,'--','Color','k');
         xline(300,'--','Color','k');
         xticks([0:200:800])
         xticklabels([-200:200:600])
         ylabel('Fr(Hz)')
         box off
         hold on;
         ax = gca;
         ax.XAxis.FontSize = 20;
         ax.YAxis.FontSize = 20;
         movegui('east');

         subplot(2,1,2)

     hit_idx = idx_k_ge{2,g}{c}(13,:)  ;
        fa_idx =  idx_k_ge{2,g}{c} (14,:) ;
        miss_idx =  idx_k_ge{2,g}{c} (15,:) ;
        cr_idx  =  idx_k_ge{2,g}{c} (16,:) ;

          hit_idx(hit_idx == 0) =nan ;
         hit_idx =  hit_idx(~isnan(hit_idx)) ;

          fa_idx(fa_idx == 0) =nan ;
         fa_idx =  fa_idx(~isnan(fa_idx)) ;

          miss_idx(miss_idx == 0) =nan ;
         miss_idx =  miss_idx(~isnan(miss_idx)) ;

          cr_idx(cr_idx == 0) =nan ;
         cr_idx =  cr_idx(~isnan(cr_idx)) ;

         color_mean = {[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};

         size_tr = [length(hit_idx),length(miss_idx), length(fa_idx),length(cr_idx)] ;
         max_size = max(size_tr) ;
        trial_types_raster = cat(1,T_k_ge{2,g}{c}(hit_idx,:),T_k_ge{2,g}{c}(miss_idx,:), T_k_ge{2,g}{c}(fa_idx,:),T_k_ge{2,g}{c}(cr_idx,:)) ;

        figure(c+1) 
     imagesc(mat2gray(trial_types_raster))

     hold on
     rectangle( 'Position' , [0 0 50 size_tr(1,1)+1],'Facecolor',color_mean{1},'EdgeColor','none')
     hold on
     yline(size_tr(1,1)+1,'--w')
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+1) 50 size_tr(1,2)+1],'Facecolor',color_mean{3},'EdgeColor','none')
     hold on
     yline((size_tr(1,1)+size_tr(1,2)),'--w')
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+size_tr(1,2)+1) 50 size_tr(1,3)],'Facecolor',color_mean{2},'EdgeColor','none')
     hold on
     yline((size_tr(1,1)+size_tr(1,2)+size_tr(1,3)+1),'--w')
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+size_tr(1,2)+ size_tr(1,3)+1) 50 size_tr(1,4)],'Facecolor',color_mean{4},'EdgeColor','none')
     hold on
xline(200,'--w');
         xline(300,'--w');
     xticks([0:200:800])
     xticklabels([-200:200:600])
     ylim([0.5 sum(size_tr)])
     yticks([1 sum(size_tr)/2 sum(size_tr) ])
     xlabel('time(ms)')
     ylabel('trial number')
% 
%      c = colorbar 
% c.FontSize = 20 
% ax = gca;
% ax.XAxis.FontSize =20;
% ax.YAxis.FontSize = 20;
% movegui('east');
% box off;
% ylabel(c,'normalized FR (Hz)')


rec_idx_th{2,g}(c)
mouse_idx_th{2,g}(c)
area_idx_th{2,g}(c)

     end
 end
  filename = 'Fig4e'
filepath = fullfile(directory, filename);
recent_figure = figure(c);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig4d'
filepath = fullfile(directory, filename);
recent_figure = figure(c+1);
saveas(recent_figure, filepath, 'svg');

    %% Figure 4 G / plot stimulus related & choice related AUC mean 
    clc
close all 
a = [1 2 3] ;
ac = a ;  % color code
% plot mean (optional include shuffled mean)
for e = 1:length(a) 
    figure(e)
    for   g = 1:numel(GNG_rec_all_cell_exp) % run per group

         length_bins = size(  GNG_analysis.mean_vec{a(e),g},2) ;
        plot([1, 1:length_bins-1],  GNG_analysis.mean_vec{a(e),g},'Color',color_eh{ac(e)},...
            'linestyle',L{g},'linewidth',2);
        hold on
        
        patch([[1,1:length_bins-1] flip([1,1:length_bins-1])] , [  GNG_analysis.mean_vec{a(e),g} + ....
              GNG_analysis.stde_vec{a(e),g} flip(  GNG_analysis.mean_vec{a(e),g} -  GNG_analysis.stde_vec{a(e),g})], color_eh{ac(e)} ,...
            'facealpha' , 0.2, 'EdgeColor','none')
        hold on
        
        plot([1, 1:length_bins-1],   GNG_analysis.mean_vec_shuf{a(e),g}+ ....
             GNG_analysis.std_vec_shuf{a(e),g},'Color',[0.5 0.5 0.5],...
            'linestyle',L{g},'linewidth',2);
        hold on
        

        xline(200,'--','Color','k');
        xline(300,'--','Color','k');

          xticks([0 6 12 18])
             xticklabels([ -200:200:400])
        xlim([1 18])
        ylim([0.5 0.6])
        xlabel('time (ms)')
        ylabel(' Discrimination (abs. AUC)')
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
        box off;
    end
    
end
  filename = 'Fig4f_easy'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
 filename = 'Fig4f_hard'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
 filename = 'Fig4f_choice'
filepath = fullfile(directory, filename);
recent_figure = figure(3);
saveas(recent_figure, filepath, 'svg');
%% Figure 4H violinplot AUC onset latency
close all
clc
clear mouse rec var
figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.onset_latency_AUC{e,1},  GNG_analysis.onset_latency_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ;
    ylim ([0 1000])
    ylabel('time(ms)')

    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
end

%% Figure 4I violinplot AUC width 

figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.width_AUC{e,1},  GNG_analysis.width_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ; 
    ylim ([0 1000])
    ylabel('time(ms)')
    
    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;

end
 
%% maximal discrimination (AUC)

figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.max_AUC{e,1},  GNG_analysis.max_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ; 
    ylim ([0.5 1.1])
    ylabel('max discrimination (AUC)')
    
    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
end 

%% LME statistics per parameter and difficulty
clc
var_onset = GNG_analysis.onset_latency_AUC ;
var_width = GNG_analysis.width_AUC ;
var_max = GNG_analysis.max_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;
for e = 1:3
    for g = 1:numel(GNG_rec_all_cell_exp)
        var_onset{e,g} = var_onset{e,g}(~isnan(var_onset{e,g})) ;
        var_width{e,g} = var_width{e,g}(~isnan(var_width{e,g})) ;
        var_max{e,g} = var_max{e,g}(~isnan(var_max{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
      
    end
end 



for e = 1:3
% Define the table with all the necessary variables
LME_AUC = table;
LME_AUC.onset = [ var_onset{e,1}'  , var_onset{e,2}']' ;
LME_AUC.width =  [ var_width{e,1}'  , var_width{e,2}']' ;
LME_AUC.maxi =    [ var_max{e,1}'  , var_max{e,2}']' ;
LME_AUC.mouse = [  mouse{e,1}'  , mouse{e,2}' + max(mouse{e,1}')]' ;
LME_AUC.rec = [ rec{e,1}'  , rec{e,2}']' ;
LME_AUC.age = [ ones(size(var_max{e,1}))' ,  ones(size(var_max{e,2}))'+ 1]' ;

% Convert variables to categorical
LME_AUC.mouse = categorical(LME_AUC.mouse);
LME_AUC.rec = categorical(LME_AUC.rec);
LME_AUC.age = categorical(LME_AUC.age);

% Fit the linear mixed effects model
formula_onset = 'onset ~ age  + (1|mouse) + (1|rec)' ;
formula_width = 'width ~ age + (1|mouse) + (1|rec)' ;
formula_max = 'maxi ~ age + (1|mouse) + (1|rec)' ;

lme_onset{e} = fitlme(LME_AUC, formula_onset);
lme_width{e}= fitlme(LME_AUC, formula_width);
lme_max{e} = fitlme(LME_AUC, formula_max);
 
disp(lme_onset{e})
disp(lme_width{e})
disp(lme_max{e})

% Extracting necessary information from the LinearMixedModel object
fixedEffects{e} = lme_onset{e}.Coefficients ;
fixedEffects{e} = [fixedEffects{e} ; lme_width{e}.Coefficients] ;
fixedEffects{e} = [fixedEffects{e} ; lme_max{e}.Coefficients] ;

randomEffects{e} = [lme_onset{e}.randomEffects lme_width{e}.randomEffects lme_max{e}.randomEffects] ;

%post_hoc comparison
fixedEffects{e}.pValue = fixedEffects{e}.pValue*3
fixedEffects{e}.pValue(fixedEffects{e}.pValue >1) = 0.9999 ;
fixedEffects{e}.Estimate = round(fixedEffects{e}.Estimate,4)
fixedEffects{e}.SE = round(fixedEffects{e}.SE,4)
fixedEffects{e}.tStat = round(fixedEffects{e}.tStat,4)
%fixedEffects{e}.pValue = round(fixedEffects{e}.pValue,4)
fixedEffects{e}.Lower = round(fixedEffects{e}.Lower,4)
fixedEffects{e}.Upper = round(fixedEffects{e}.Upper,4)

% Create a table for fixed effects
fixedEffectsTable = table(fixedEffects{e}.Estimate, fixedEffects{e}.SE, fixedEffects{e}.tStat, fixedEffects{e}.DF, fixedEffects{e}.pValue,  fixedEffects{e}.Lower,...
    fixedEffects{e}.Upper,...
    'VariableNames', {'Estimate', 'StandardError', 'tStatistic', 'DF', 'pValue', 'CI lower', 'CI upper'});

writetable(fixedEffectsTable,'LME_AUC_fixed.xlsx','Sheet',e)

end 
%% plot 3d per stimulus or chocie 
close all
for g = 1:2
    for e = 1:3
        figure
        subplot(1,2,e)
        plot3(GNG_analysis.onset_latency_AUC{e,g},GNG_analysis.width_AUC{e,g},GNG_analysis.max_AUC{e,g}...
            ,'LineStyle','none','Marker','o','Color',color_eh{e},'MarkerFaceColor',color_eh{e})
        hold on
        grid on
        xlabel('onset')
        xlim([0 1000])
        ylabel('width')
        ylim([0 1000])
        zlabel('max')
        zlim([0.5 1])
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        ax.ZAxis.FontSize = 20;
        movegui('east');
        box off;
    end
end

  filename = 'Fig4g_easy'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig4g_hard'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig4g_choice'
filepath = fullfile(directory, filename);
recent_figure = figure(3);
saveas(recent_figure, filepath, 'svg');

%% supplementary figure 53 plot 3d per area 
close all
Colors_area = {[.1 .3 .8], [.5 .4 .9], [0 .5 .6],[0.4940 0.1840 0.5560]};
for g = 1:2
    figure
    for e = 1:3
        idx_AUDd = (find(area_idx_th{e,1} == 1 ))
        idx_AUDp = (find(area_idx_th{e,1} == 2 ))
        idx_AUDv = (find(area_idx_th{e,1} == 3 ))
        idx_TEa = (find(area_idx_th{e,1} == 4 ))

        subplot(1,3,e)
        plot3(GNG_analysis.onset_latency_AUC{e,g}(idx_AUDd),GNG_analysis.width_AUC{e,g}(idx_AUDd),GNG_analysis.max_AUC{e,g}(idx_AUDd)...
            ,'LineStyle','none','Marker','o','Color',Colors_area{1},'MarkerFaceColor',Colors_area{1})
        hold on
        plot3(GNG_analysis.onset_latency_AUC{e,g}(idx_AUDp),GNG_analysis.width_AUC{e,g}(idx_AUDp),GNG_analysis.max_AUC{e,g}(idx_AUDp)...
            ,'LineStyle','none','Marker','o','Color',Colors_area{2},'MarkerFaceColor',Colors_area{2})
        hold on
        plot3(GNG_analysis.onset_latency_AUC{e,g}(idx_AUDv),GNG_analysis.width_AUC{e,g}(idx_AUDv),GNG_analysis.max_AUC{e,g}(idx_AUDv)...
            ,'LineStyle','none','Marker','o','Color',Colors_area{3},'MarkerFaceColor',Colors_area{3})
        hold on
        plot3(GNG_analysis.onset_latency_AUC{e,g}(idx_TEa),GNG_analysis.width_AUC{e,g}(idx_TEa),GNG_analysis.max_AUC{e,g}(idx_TEa)...
            ,'LineStyle','none','Marker','o','Color',Colors_area{4},'MarkerFaceColor',Colors_area{4})

        grid on
        xlabel('onset')
        xlim([0 1000])
        ylabel('width')
        ylim([0 1000])
        zlabel('max')
        zlim([0.5 1])
    end
end 
  filename = 'FigS43_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
  filename = 'FigS43_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
