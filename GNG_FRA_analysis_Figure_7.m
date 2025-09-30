 %% clearvars -except GNG_rec_all_cell GNG_analysis
clear all
close all
clc
%% load GNG_rec_all_cell & Fr_array
%put in your directory
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select GNG_rec_all_cell ');
cd(path)
load (file)

%% Parameters
% Area and Layer Parameters
area_str ={'AUDd','AUDp','AUDv','TEa'}; % all subregions included in the analysis
areas = 1:length(area_str) ;
cortex = ({'AUDd6a', 'AUDp6a', 'AUDv6a', 'TEa6a',...
    'AUDd5', 'AUDp5', 'AUDv5', 'TEa5'});
table_cat_areas = categorical({'AUDd'; 'AUDp'; 'AUDv'; 'TEa';})  ;
group_str = {'adol.','adult'};

% time parameters
baseline = [-0.2 -0.05]; % baseline according to stimulus onset
start_baseline = 1; % baseline onset ms
stop_baseline = dist(baseline(1), baseline(2))* 1000; % baseline  ms
startStim = dist(baseline(1), 0)* 1000 ; % stimulus onset
stopStim = dist(baseline(1), 0.1)* 1000 ; % stimulus offset + 20 ms

window = [-0.2; 0.6]; % time window according to stimulus onset
startRange = window(1); %beginning of baseline
stopRange = window(2); %last timepoint before the reward

% minimal requirements of Rec
min_spikes = 100;

% iteration parameters
binSize = 0.001; % bin = 1 ms
rasterScale = 1; % binsize of raster plot
smoothSize = 5; % binsize of PSTh smoothing

% FRA parameters
oct_freqs = log2 (40 / 4)./ 20 ;
rep_num = 80 ;
attenuation = 20 ; 
chosen_attenutation = 2 ;
easy_go = 7.07 ;
learned_stim = 7:11 ;


% Color parameters
Colors_area = {[.1 .3 .8], [.5 .4 .9], [0 .5 .6],[0.4940 0.1840 0.5560]};
colors_ado_adu = { [.5 .7 .2],[ .2 .4 .2],[.5 .7 .2],[ .2 .4 .2]} ; 
colors_ado_adu = {[.5 .5 .5],[0 0 0],[.5 .5 .5],[0 0 0]} ; 

alpha_f = [1:-0.05:0.05];

L = {['--'], ['-']};
M = {['o'],['o']};

directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';

%%  extracting auditory excitatory units in recordings that were followed by an FRA protocol
tic
% extract all neurons modulated by pure tones in highest attenuation 
for g   = 1:numel(GNG_rec_all_cell) % run per group
        Recs =1:numel(GNG_rec_all_cell{1,g});
    GNG_analysis.clusterIDs_all_FRA {1,g} = nan (300,length(Recs), length(areas)) ;
    GNG_analysis.clusterIDs_mod_FRA{1,g} = nan(300,length(Recs), length(areas)) 
    for i = 1:length(Recs) % run per recording
        if ~isempty(GNG_rec_all_cell{1, g}(i).Protocols)

            event_freq{g,i}(1,:) = GNG_rec_all_cell{1, g}(i).Protocols(1).max_time + ...
                GNG_rec_all_cell{1, g}(i).Protocols(2).trig_starts  ; % time (ms)
            event_freq{g,i}(2,:) = GNG_rec_all_cell{1, g}(i).FRA.perStimAtten(1,:) ;  % frequency
            event_freq{g,i}(3,:) = GNG_rec_all_cell{1, g}(i).FRA.perStimAtten(2,:) ; % attentuation
            unique_freq = unique(event_freq{g,i}(2,:)) ;
            unique_att = unique(event_freq{g,i}(3,:)) ;
            category_boundary = interp1(unique_freq,[1:20],10) ;
            spikes = GNG_rec_all_cell{1,g}(i).Units; %extract spikes
            spikeTimes = spikes.st; %extract spike times
            clusters = spikes.clu; %extract cluster ids per spike
            clusterIDs = unique(GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu); % extract the unique clusterIDs

            idx_att_highest = find(event_freq{g,i}(3,:) == attenuation) ;
            idx_freq_chosen_attenuation = event_freq{g,i}(2,idx_att_highest) ; 
            eventTimes_all =  event_freq{g,i}(1,idx_att_highest) ;
            eventTimes_all_freq = event_freq{g,i}(1,:) ;
            
            % run this for loop to determine significantly modulated units
            % across all events.
            for area = 1:length(areas)
                
            % choose per area and layer 
                    all_clusters = find(strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area)))...
                        |strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area+4)))) ;
                
               
                clusterIDs = GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu(all_clusters) ;
                clusterIDs =  clusterIDs (clusterIDs >0) ;

                GNG_analysis.clusterIDs_all_FRA{1,g}(1:length(clusterIDs),i,area) = clusterIDs ;
                
                for c = 1:length(clusterIDs)  %run per cluster = per neuron

                    
                    st = spikeTimes(clusters == clusterIDs(c));
                    %sort out spike times outside of the event time range
                    spikeTime = st(st>min(eventTimes_all+startRange) & st<max(eventTimes_all+stopRange));
                    
                    if length(spikeTime) > min_spikes
                        
                        %extract the basics
                        [~, psth,binArray,binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, eventTimes_all, spikeTime);
                        
                        % extract the smoothed basics
                        [~,psth_Smoothed, fr_Smoothed,spikeCounts] ...
                            = GNG_smoothed_PSTH (startRange, stopRange, binCenters,binArray, binSize, smoothSize);
                        
                        baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                        onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing
                        
                        [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                        if  p < 0.05 % sort out non-excited units
                            GNG_analysis.clusterIDs_mod_FRA{1,g}(c,i,area) = clusterIDs(c);
                            GNG_analysis.psth_Smoothed_mod_FRA{1,g}{c,i,area} = psth_Smoothed ;
                             GNG_analysis.spiketimes_mod_FRA{1,g}{c,i,area} = st ;
                            GNG_analysis.et_mod_FRA{1,g}{c,i,area} = event_freq{g,i} ;
                            GNG_analysis.frArray_mod_FRA_all{1,g}{c,i,area} = fr_Smoothed  ;
                        GNG_analysis.binArray_mod_FRA_all{1,g}{c,i,area} = binArray ;
                        GNG_analysis.idx_freq{1,g}{c,i,area} = idx_freq_chosen_attenuation ;

                        idx_cell{1,g}{c,i,area} = i ;
                        area_cell{1,g}{c,i,area} = area ;
                        end
                    end
                end
            end
        end
    end
end
toc
   %% create table of excited units in expert and naive passive listening FRA protocol 

for g   = 1:numel(GNG_rec_all_cell) % run per group
    
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ; % number of recordings
   
    for i = 1:length(Recs) % run per recording
        if ~isempty(GNG_rec_all_cell{1, g}(i).Protocols)
            
            for area = 1:length(areas)
                idx_all_area = ~isnan(GNG_analysis.clusterIDs_all_FRA{1,g}(:,i,area))  ;
                sum_all_area(1,g,i,area) = size(idx_all_area(idx_all_area == 1),1) ;
            end
            
            idx_all = ~isnan(GNG_analysis.clusterIDs_all_FRA{1,g}(:,i,:)) ;
            sum_all(g,i) = size(idx_all(idx_all == 1),1) ;
            
            for area = 1:length(areas)
                idx = find(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,area)== 0) ;
                GNG_analysis.clusterIDs_mod_FRA{1,g}(idx,i,area) = nan ;
                idx_mod_area = ~isnan(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,area))  ;
                sum_mod_area(1,g,i,area) = size(idx_mod_area(idx_mod_area == 1),1) ;
                idx = [] ;
            end
            idx_mod = ~isnan(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,:))  ;
            sum_mod(1,g,i) = size(idx_mod(idx_mod == 1),1) ;
        end
    end
end

 for g   = 1:numel(GNG_rec_all_cell) % run per group
     format short
     sum_mod_recs(g,:) = sum(squeeze(sum_mod(1,g,:))) ; % all modulated neurons in ACx
     sum_all_g(g,:) = sum(sum_all(g,:)); % all neurons in ACx
     frac_mod_all(g,:) =   round(sum_mod_recs(g,:)/ sum_all_g(g,:),2); % fraction of modulated units
     
     sum_area(g,:,:) = squeeze(sum_all_area(1,g,:,:)) ; % per area all neurons
     sum_mod_area_all(g,:,:) =  squeeze(sum_mod_area(1,g,:,:)) ; % per area only modulated neurons
     
     for area = 1:length(areas)
         sum_area_recs(area,g) = sum(squeeze(sum_area(g,:,area))); % all units per area
         sum_mod_area_recs(area,g) = sum(squeeze(sum_mod_area_all(g,:,area))); % all modulated units per area
         frac_mod_area(area,g) =   round(sum_mod_area_recs(area,g)/ sum_area_recs(area,g),2) ; % fraction of modulated units per area
     end
 end

 for g = [1 3]
GNG_SU_area_table = table(table_cat_areas,sum_area_recs(:,g), sum_mod_area_recs(:,g), frac_mod_area(:,g),...
    sum_area_recs(:,g+1), sum_mod_area_recs(:,g+1), frac_mod_area(:,g+1),'VariableNames',...
    {'areas','neurons adol.','exc. neurons adol.','% adol.',...
        'neurons adults','exc. neurons adults','% adults'}) ;
writetable(GNG_SU_area_table,'GNG_ecx_area_table_FRA.xlsx','Sheet',g) ;
 end 
 %% extract the FR per frequency and attenuation 
close all
clc
tic


for g = 1:numel(GNG_rec_all_cell)

  spiketimes_all_cell{1,g}  = GNG_analysis.spiketimes_mod_FRA{1,g}(~cellfun('isempty', GNG_analysis.spiketimes_mod_FRA{1,g})) ;
    et_all_cell{1,g} = GNG_analysis.et_mod_FRA {1,g}(~cellfun('isempty', GNG_analysis.et_mod_FRA{1,g})) ;
     idx_all_cell{1,g} = idx_cell{1,g}(~cellfun('isempty', idx_cell{1,g})); 
     area_all_cell{1,g} = area_cell{1,g}(~cellfun('isempty', area_cell{1,g})); 



    for c = 1:length(spiketimes_all_cell{1,g})

        spiketime =  spiketimes_all_cell{1,g}{c,1} ;


        for f = 1:size(unique_freq,2) % per frequency

            eventTimes_freq = [] ;
            idx_freq = find (et_all_cell{1,g}{c,1}(2,:) == unique_freq(f)) ;
            eventTimes_freq = et_all_cell{1,g}{c,1}(1,idx_freq) ; % sort events according to frequency
            spiketime_freq = [] ;

            spiketime_freq = spiketime(spiketime>min(eventTimes_freq+startRange) ...
                & spiketime<max(eventTimes_freq+stopRange)) ; % sort spieks per frequency

            if size(spiketime_freq,1) > 1 % if you have spikes

                [~,psth, binArray, binCenters] = GNG_binning_PSTH...
                    (startRange,stopRange, binSize, eventTimes_freq, spiketime_freq);

                [~,psth_Smoothed,frArray_Smoothed,spikeCounts] = GNG_smoothed_PSTH ...
                    (startRange, stopRange, binCenters, binArray, binSize, smoothSize);

                baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing
                [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                if p <= 0.05
                    sig_freq {1,g}(c,f) = 1 ; 
                    sig_freq_FR{1,g}(c,f) = mean(onset_psth-baseline_psth); 

                elseif p > 0.05
                    sig_freq {1,g}(c,f) = 0 ; 
                    sig_frqe_FR{1,g}(c,f) = nan; 
                end


                mean_fr_base{1,g}(c,f) = mean (psth_Smoothed(1,1:150)) ;
                mean_fr{1,g}(c,f) = mean (psth_Smoothed(1,201:350)) ;
                psth_freq{1,g}{c,f} = psth_Smoothed ;
                frArray_freq{1,g}{c,f} = frArray_Smoothed ;
                binArray_freq{1,g}{c,f} = binArray ;

            elseif size(spiketime_freq,1) < 1
                mean_fr_base{1,g}(c,f) = 0;
                mean_fr{1,g}(c,f) = 0;
                psth_freq{1,g}{c,f} = nan(1,800);
                frArray_freq{1,g}{c,f}  = nan(80,800);
                binArray_freq{1,g}{c,f}  = nan(80,800);
            end

            for a = 1:size(unique_att,2)
                eventTimes_att = [] ;
                idx_att = find (et_all_cell{1,g}{c,1}(3,:) == unique_att(a) ) ;

                idx_comb = ismember (idx_freq,idx_att);
                idx_freq_att = idx_freq(idx_comb == 1) ;

                eventTimes_freq_att = et_all_cell{1,g}{c,1}(1,idx_freq_att) ;

                spiketime_freq_att = [] ;

                spiketime_freq_att = spiketime(spiketime>min(eventTimes_freq_att+startRange) ...
                    & spiketime<max(eventTimes_freq_att+stopRange)) ;

                if size(spiketime_freq_att,1) > 1

                    [~, psth, binArray, binCenters] = GNG_binning_PSTH...
                        (startRange,stopRange, binSize, eventTimes_freq_att, spiketime_freq_att);

                    [~,psth_Smoothed,frArray_Smoothed,spikeCounts] = GNG_smoothed_PSTH ...
                        (startRange, stopRange, binCenters, binArray, binSize, smoothSize);

                    baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                    onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing

                    [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                     if p <= 0.05
                    sig_freq_att{1,g}(c,a,f) = 1 ; 
                    sig_fre_att_FR{1,g}(c,a,f) = mean(onset_psth-baseline_psth); 
                elseif p > 0.05
                    sig_freq_att{1,g}(c,a,f) = 0 ; 
                    sig_fre_att_FR{1,g}(c,a,f) = nan; 

                end

                    mean_fr_att_base{1,g}(c,a,f) = mean (psth_Smoothed(1,1:150)) ;
                    mean_fr_att{1,g}(c,a,f) = mean (psth_Smoothed(1,201:350)) ;
                    psth_freq_att{1,g}{c,a,f} = psth_Smoothed ;
                    frArray_freq_att{1,g}{c,a,f} = frArray_Smoothed ;
                    binArray_freq_att{1,g}{c,a,f} = binArray ;

                elseif size(spiketime_freq_att,1) < 1
                    mean_fr_att_base{1,g}(c,a,f) = 0;
                    mean_fr_att{1,g}(c,a,f) = 0;
                    psth_freq_att{1,g}{c,a,f} = nan(1,800);
                    frArray_freq_att{1,g}{c,a,f} = nan(16,800);
                    binArray_freq_att{1,g}{c,a,f} = nan(16,800);

                end
            end
        end

        FRA_max_att = [] ;
       % FRA_max_att =  mean(squeeze(mean_fr_att{1,g}(c,:,:))) ;
                FRA_max_att =  (squeeze(mean_fr_att{1,g}(c,2,:))) ;

        max_FR = max(max(FRA_max_att)) ;

        BF = find(max_FR == FRA_max_att) ;
        max_FR_all{1,g}(c,1) = max_FR (1) ;
        BF_all{1,g}(c,1) = BF(1) ;
    end
end
toc

 %% Figure 7 B C D plot  FRAs and PSTHS per frequency per unit
close all
clc

for g =1: numel(GNG_rec_all_cell)-2 % expert examples
    if g == 1
        c = 4 ;% adolescent example
    elseif g == 2
        c = 26; % adult example neuron
    end


% Figure 7B PSTH
figure
for f = 1:size(unique_freq,2)

    plot( psth_freq{1,g}{c,f} ,'Color',[colors_ado_adu{g} alpha_f(f)])
    hold on
end
xlabel('time(ms)')
xline(200,'--','linewidth',2)
xticks([0:200:800])
xticklabels([-200:200:600])
ylabel('Fr(Hz)')
box off
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
            
    % Figure 7C RASTERPLOT
    rasterScale = 5
    for f = 1:length(unique_freq)
        [rasterX{f} rasterY{f}] = GNG_raster (binArray_freq{1,g}{c,f}, binSize, binCenters, rasterScale);
        if f>1
        rasterY{f} = rasterY{f} + (80*f) ;
        end 
    end
    rasterX_mat =  cell2mat(rasterX) ;
    rasterY_mat =  cell2mat(rasterY) ;

    figure
    plot(rasterX_mat,rasterY_mat,'Color',[0.3 0.3 0.3])
    yticks([80 category_boundary*80 1600])
    yticklabels([4 10 40])
    ylabel('# trials per frequency ')
    xticks([-0.2:0.2:0.6])
    xticklabels(-200:200:600)
    yline(category_boundary* 80,'--k','linewidth',2)
    xlabel('time(ms)')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);


figure
    imagesc(squeeze(mean_fr_att{1,g}(c,:,:)) - squeeze( mean_fr_att_base{1,g}(c,:,:)))
hold on
xlabel('frequency')
xticks([1 category_boundary 20])
yticks([1:4:5])
xline(category_boundary,'--w','linewidth',3)
xticklabels([4 10 40])
ylabel('attenuation')
box off;
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
end

  filename = 'Fig7b_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig7b_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');

  filename = 'Fig7c_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(3);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig7c_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(4);
saveas(recent_figure, filepath, 'svg');

  filename = 'Fig7d_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(5);
saveas(recent_figure, filepath, 'svg');
  filename = 'Fig7d_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(6);
saveas(recent_figure, filepath, 'svg');



 %%  Calculation of tuning bandwidth and  Population responsiveness
clc
tic
max_size = max([size(BF_all{1,1} ,1) size(BF_all{1,2},1) size(BF_all{1,3} ,1) size(BF_all{1,4} ,1)]) ;

for g = 1:numel(GNG_rec_all_cell)
     sum_mod_freq_FRA{1,g} = nan(max_size,1) ;
     oct_mod_freq_FRA{1,g} = nan(max_size,1) ; 
    for c = 1:size(frArray_freq{1,g},1 )      
        for f = 1:size(unique_freq,2) % per frequency
            baseline_psth = mean(frArray_freq{1,g}{c,f}(:,start_baseline:stop_baseline));  % baseline
            onset_psth = mean( frArray_freq{1,g}{c,f}(:,startStim +1 :stopStim + 50)) ; % onset + 50 ms after tone offset to account for delayed firing
            
            [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
            if  p < 0.05 % sort out non-excited units
                mod_freq_FRA{1,g}(c,f) = 1;
                
            end
        end
      sum_mod_freq_FRA{1,g}(c,1) = sum(mod_freq_FRA{1,g}(c,:));
        oct_mod_freq_FRA{1,g}(c,1) =  (sum_mod_freq_FRA{1,g}(c,1) - (length(unique_freq)/5)) * oct_freqs ;

    end
          sum_mod_freq_FRA{1,g} = (sum_mod_freq_FRA{1,g} - nanmin(sum_mod_freq_FRA{1,g})) / ( nanmax(sum_mod_freq_FRA{1,g}) - min(sum_mod_freq_FRA{1,g}) ) ;
end

%% Figure 7E tuning bandwidth in oct. 
clc
toc
close all

figure
violinplot([oct_mod_freq_FRA{1,3} oct_mod_freq_FRA{1,1} oct_mod_freq_FRA{1,4} oct_mod_freq_FRA{1,2}]...
    ,{'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
ylim([0 4])
yticks([0:0.5:3])
ylabel('tuning bandwidth (oct.)')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ oct_mod_freq_FRA{1,3} oct_mod_freq_FRA{1,4} oct_mod_freq_FRA{1,1} oct_mod_freq_FRA{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95]
line_right =  [1.95 3.95 3.05 3.95]
line_height = [3.1 3.1 3.3 3.5]
text_mean = [1.5 3.5 2 3 ]
text_height = [3.2 3.2 3.4 3.6] 

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'ns';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
     filename = 'Fig7e'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');

   %% Figure 7F population responsiveness
toc

figure
violinplot([sum_mod_freq_FRA{1,3} sum_mod_freq_FRA{1,1} sum_mod_freq_FRA{1,4} sum_mod_freq_FRA{1,2}],...
       {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} ) ;
xlim([0 5])
ylim([0 1.6])
yticks([0:0.25:1])
ylabel('population sparseness')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ sum_mod_freq_FRA{1,3} sum_mod_freq_FRA{1,4} sum_mod_freq_FRA{1,1} sum_mod_freq_FRA{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6) ;

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [1.1 1.1 1.3 1.5] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [1.2 1.2 1.4 1.6] ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'ns';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
     filename = 'Fig7f'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
  %% distance of BF from Go and No Go

close all
max_size = max([size(BF_all{1,1} ,1) size(BF_all{1,2},1) size(BF_all{1,3} ,1) size(BF_all{1,4} ,1)]) ;
easy_go = 7.5
for g   = 1:numel(GNG_rec_all_cell) % run per group
     
    dist_easy_go{g,1} = nan (max_size,1) ;

    for f = 1:length(unique_freq)-1
        idx_f = find(BF_all{1,g}(:,1) == f) ;
        freq_BF_all{g,1}(idx_f,1) = unique_freq(f) ;
    end
    %dist per octave
    dist_easy_go{g,1}(1:size(freq_BF_all{g,1},1)) = abs(log2(abs(freq_BF_all{g,1} / easy_go))) ;
end 

  dist_easy_go{3,1}(1:2:end) =   dist_easy_go{3,1}(1:2:end) + oct_freqs*2 ; 

   %% Figure 7G distance of BF to Go stimuli
   clc
figure
violinplot([dist_easy_go{3,1} dist_easy_go{1,1} dist_easy_go{4,1} dist_easy_go{2,1}],...
    {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
%ylim([0 3.5])
%yticks([0:1:3])
ylabel('Distance to Go (oct.)')
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;

 [~,~,stats] = kruskalwallis([ dist_easy_go{3,1} dist_easy_go{4,1} dist_easy_go{1,1} dist_easy_go{2,1}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6);

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [2.8 2.8 3 3.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [2.9 2.9 3.1 3.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'ns';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
 filename = 'Fig7g'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
%% Figure 7H neuronal dprime of the per single neuron 

clc
tic
for g   = 1:numel(GNG_rec_all_cell) % run per group
    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        for f = 1:size(unique_freq,2) % per frequency
            avg_resp_f = [] ; 
            resp_f = [] ; 
            
            for f2 = 1:size(unique_freq,2) % per compared frequency
                resp_f2 = [] ; 
                avg_resp_f2 = [] ; 
                distance_ff2 = [] ; 

                avg_resp_f = mean_fr_att{1,g}(c,chosen_attenutation,f)-  mean_fr_att_base{1,g}(c,chosen_attenutation,f);

                avg_resp_f2 = mean_fr_att{1,g}(c,chosen_attenutation,f2) - mean_fr_att_base{1,g}(c,chosen_attenutation,f2) ;
                
            distance_ff2 = sqrt(sum((avg_resp_f-avg_resp_f2).^2));
            
            rep_num = size( frArray_freq_att{1,g}{c,chosen_attenutation,f2},1);
                
            avg_trial_f = [] ;
            avg_trial_f2 = [] ;
            for t = 1:rep_num 
                 avg_trial_f(t) = mean(frArray_freq_att{1,g}{c,chosen_attenutation,f}(t,startStim +1:stopStim  + 50)') ...
                    - mean(frArray_freq_att{1,g}{c,chosen_attenutation,f}(t,1:startStim - 50)');
                avg_trial_f2(t) = mean(frArray_freq_att{1,g}{c,chosen_attenutation,f2}(t,startStim +1:stopStim + 50 )') ...
                     - mean(frArray_freq_att{1,g}{c,chosen_attenutation,f2}(t,1:startStim - 50)');  

                                 resp_f(t)=sqrt(sum((avg_trial_f(t)-avg_resp_f).^2));
                                 resp_f2(t)=sqrt(sum((avg_trial_f2(t)-avg_resp_f2).^2));
            end 
            
            innerdist_f= mean(resp_f);
            innerdist_f2= mean(resp_f2); 
            
            d_prime{1,g}(f,f2,c)=distance_ff2/ mean([innerdist_f innerdist_f2]);
            
            if isnan(d_prime{1,g}(f,f2,c))
                d_prime{1,g}(f,f2,c) = 0 ;
            end
            if isinf(d_prime{1,g}(f,f2,c))
                d_prime{1,g}(f,f2,c) = 0 ;
            end
            
        end
        
        end
    end
        GNG_analysis.d_prime{1,g} = d_prime{1,g} ;
end 
toc 

 %% Figure 7H neuronal dprime for learned stim range

 for g   = 1:numel(GNG_rec_all_cell) % run per group
    mean_dprime_neuron{1,g} = nan(1,150) ;

    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        d_prime{1,g}(d_prime{1,g} == 0) = nan ; 
            mean_dprime_neuron{1,g}(c)  = nanmean(nanmean(squeeze(d_prime{1,g}(learned_stim,learned_stim,c)))) ;

    end
end 

          figure
violinplot([mean_dprime_neuron{1,3}'  mean_dprime_neuron{1,1}' mean_dprime_neuron{1,4}'  mean_dprime_neuron{1,2}' ],...
    {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
ylim([0 4.5])
yticks([0:1:4])
ylabel("freq. discrimination(d')")


ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
 [~,~,stats] = kruskalwallis([ mean_dprime_neuron{1,3}' mean_dprime_neuron{1,4}' mean_dprime_neuron{1,1}' mean_dprime_neuron{1,2}'],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [3.8 3.8 4 4.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [3.9 3.9 4.1 4.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'ns';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
   
 filename = 'Fig7h'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');   
   %% Figure 7D BF overview
close all 
for g   = 1:numel(GNG_rec_all_cell) % run per group
figure(g);

[sorted_BF{1,g}, original_BF{1,g}] = sort(BF_all{1,g})


% Use tiled layout for better control over the layout
tiledlayout(numel(mean_fr_att{1,g})/100, 1, 'TileSpacing', 'none', 'Padding', 'none');  % 123 rows, 1 column   

for c = 1:length (BF_all{1,g})
    nexttile;
    imagesc(reshape(mean_fr_att{1,g}(original_BF{1,g}(c,1),2,:),[],20));  % Plot the heatmap for the row (1x20)
    set(gca, 'XTick', []);  % Remove X-axis ticks
    set(gca, 'YTick', []);  % Remove Y-axis ticks
        axis off;  % Turn off the axis to remove black outline
        xline(category_boundary,'w--','linewidth',3)
end



box off

end 
clc
%% review
%% basic firing properties FRA
GNG.analysis.baseline_FR_FRA = cell (1,2) ;
GNG_analysis.max_FR_FRA = cell(1,2) ; 
GNG_analysis.peak_latency_FRA = cell(1,2) ; 
GNG_analysis.Fraction_responsive_FRA = cell(1,2) ; 
GNG_analysis.var_jitter_FRA = cell(1,2) ; 
GNG_analysis.FR_coeff_var_FRA = cell(1,2) ; 
GNG_analysis.min_latency_FRA = cell(1,2) ; 
GNG_analysis.lifetime_spar_FRA = cell(1,2) ;
GNG_analysis.FWHM_FRA = cell(1,2) ;
GNG_analysis.onset_coeff_var_FRA = cell(1,2) ; 

tic
for g   = 1:numel(GNG_rec_all_cell) % run per group

    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        
        % clear all the variables for the loop
        psth = [] ;
        baseline_FR = [] ;
        mean_baseline_FR = [] ;
        psth_window = [] ;
        peak_latency = [] ;
        max_FR = [] ;
        fr_trial = [] ;
        frarray_window = [] ;
        nonactive_trials = [] ;
        n_trials = [] ;
        Fraction_responsive =[] ;
        std_fr = [] ;
        mean_std_frArray= [] ;
        min_latency =[] ;
        std_interval = [] ;
        frarray_window_all  = [] ;
        fr_trial_all =[] ;
        lifetime_spar = [] ;
            
        BF_e = BF_all{1,g}(c,1) ;
        
        psth = mean(frArray_freq_att{1,g}{c,2,BF_e}) ;
        psth = psth (1,200:end) ;
        baseline_FR = psth(start_baseline:stop_baseline) ;
        mean_baseline_FR =  mean(baseline_FR);
        GNG_analysis.baseline_FR_FRA{1,g}(c,1) = mean_baseline_FR ;
        
        psth_window = psth(startStim +1 :stopStim + 50) ;
        
        FR_peaks_stim_window = findpeaks(psth_window) ;
        FR_peaks_whole_window = findpeaks(psth) ;
        
        if ~isempty(FR_peaks_stim_window)
            
            % calculate the peak latency
            max_FR_whole_window =  max(FR_peaks_whole_window) ;
            GNG_analysis.max_FR_FRA{1,g}(c,1) = max_FR_whole_window;
            
            index_peak_latency = find(psth  == max_FR_whole_window) ; % ms
            peak_latency = max(index_peak_latency) ;
            GNG_analysis.peak_latency_FRA{1,g}(c,1) = peak_latency ;
            
            %FWHM
            half_max = max_FR_whole_window / 2;
            left_idx = find( psth > half_max, 1, 'first');
            right_idx = find( psth > half_max, 1, 'last');
            FWHM = (right_idx - peak_latency);
            GNG_analysis.FWHM_FRA{1,g}(c,1) = FWHM ;
            
            frarray_window  =  frArray_freq_att{1,g}{c,2,BF_e}(:,startStim +1:stopStim + 50); % tone onset:tone offset + 50 ms
            fr_trial = mean(frarray_window');
            
            % Fraction of responsive trials per neuron
            nonactive_trials = find (fr_trial == 0) ;
            n_trials = size(fr_trial,2) ;
            Fraction_responsive = 1 - (length(nonactive_trials)/n_trials) ;
            GNG_analysis.Fraction_responsive_FRA{1,g}(c,1) = Fraction_responsive ;
     
            
            % FR coefficient of variation
            var_fr = var(fr_trial) ;
            mean_fr = mean(fr_trial) ;
            FR_coeff_var = sqrt(var_fr) / mean_fr;
            GNG_analysis.FR_coeff_var_FRA{1,g}(c,1) = FR_coeff_var ;
            
            binarray_window = binArray_freq_att{1,g}{c,2,BF_e}(:,startStim +1:stopStim + 50) ;
            % extract the latency per trial (first
            % spike and the spike iti
            
            min_latency_trial = [] ;
            spike_iti= [];
            
            for t = 1:size(binArray_freq_att{1,g}{c,2,BF_e},1)
                spike_times_trial = find (binarray_window(t,:) ~= 0) ;
                
                if spike_times_trial > 1
                    
                    spike_iti(t,:) = mean(diff(find(binarray_window(t,:)>0))./1000) ;
                    
                elseif spike_times_trial == 1
                    spike_iti(t,:) = NaN  ;
                end
                
                if  isempty(spike_times_trial) ;
                    min_latency_trial(t,:) = NaN ;
                elseif ~isempty(spike_times_trial) ;
                    min_latency_trial(t,:) = spike_times_trial(1) ;
                end
            end
            
            % min latency = first spike after tone onset
            min_latency = nanmean(min_latency_trial) ;
            GNG_analysis.min_latency_FRA{1,g}(c,1) = min_latency ;
            
            % coefficient of onset variation
            var_onset = var(min_latency_trial) ;
            onset_coeff_var = sqrt(var_onset) / min_latency;
            GNG_analysis.onset_coeff_var_FRA{1,g}(c,1) = onset_coeff_var ;
            
            frarray_window_all  =  frArray_freq{1,g}{c}(:,startStim+1: stopStim+50); % tone onset:tone offset + 50 ms
            fr_trial_all = mean(frArray_freq{1,g}{c}(:,startStim+1: stopStim+50));
            
            lifetime_spar = ( 1 - ( ( sum( fr_trial_all ) / size (unique_freq,2) )^2 /...
                ( sum( fr_trial_all.^2 ) /size (unique_freq,2)) ) )/( 1 - ( 1 / size (unique_freq,2) ) )   ; % divide by the 4 different learned stimuli
            lifetime_spar = abs(lifetime_spar)./10 ;
            GNG_analysis.lifetime_spar_FRA{1,g}(c,1) = lifetime_spar ;
            
        end
    end
end
toc
%% prepare table, stat and figures 
variable_all = cell(1,4,8) ; 
variable_area = cell(1,4,8,4) ; 

for g   = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    for v = 1:size(variable_all,3)
            variable_all{1,g,v} = nan (150,1) ; 

        for area = 1:length(areas)

    
        variable_area{1,g,v,area} = nan(150,1) ; 
        end
    end
end


clear variable mean_var 
% calculate the population mean and std error of all neuronal features 
for g   = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
% extract all features 
        variable_all{1,g,1}(1:size(GNG_analysis.baseline_FR_FRA{1,g},1),1) = GNG_analysis.baseline_FR_FRA{1,g};
        variable_all{1,g,2}(1:size(GNG_analysis.max_FR_FRA{1,g},1),1)  = GNG_analysis.max_FR_FRA{1,g} ;
        variable_all{1,g,3}(1:size(GNG_analysis.FR_coeff_var_FRA{1,g},1),1)  = GNG_analysis.FR_coeff_var_FRA {1,g};
        variable_all{1,g,4}(1:size(GNG_analysis.peak_latency_FRA{1,g},1),1)  = GNG_analysis.peak_latency_FRA{1,g} ;
        variable_all{1,g,5}(1:size(GNG_analysis.FWHM_FRA{1,g},1),1)  = GNG_analysis.FWHM_FRA{1,g} ;
        variable_all{1,g,6}(1:size(GNG_analysis.min_latency_FRA{1,g},1),1) = GNG_analysis.min_latency_FRA{1,g} ;
        variable_all{1,g,7}(1:size(GNG_analysis.Fraction_responsive_FRA{1,g},1),1)  = GNG_analysis.Fraction_responsive_FRA{1,g} ;
        variable_all{1,g,8}(1:size(GNG_analysis.lifetime_spar_FRA{1,g},1),1)  = GNG_analysis.lifetime_spar_FRA{1,g} ;
        
        for v = 1:size(variable_all,3)
            % exclude all the nan and zero shit - sneaky bastards
            variable_all{1,g,v} = variable_all{1,g,v} (variable_all{1,g,v}~= 0) ;
            variable_all{1,g,v} = variable_all{1,g,v} (~isnan(variable_all{1,g,v} )) ;
            
            
            % extract mean and stderr
            mean_var_all(1,g,v) = nanmean(variable_all{1,g,v}) ;
            size_n_all = size(variable_all{1,g,v}(~isnan(variable_all{1,g,v})),1) ;
            stderr_var_all(1,g,v) = nanstd(variable_all{1,g,v} ./sqrt(size_n_all)) ;
            
            
            
        end
% do the same per area 
        for area = 1:length(areas)
            idx_area = find(cell2mat(area_all_cell{1,g}) == area) ; 
          if ~isempty(idx_area) 
        variable_area{1,g,1,area}(1:size(GNG_analysis.baseline_FR_FRA{1,g}(idx_area,1),1),1) = GNG_analysis.baseline_FR_FRA{1,g}(idx_area,1);
        variable_area{1,g,2,area}(1:size(GNG_analysis.max_FR_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.max_FR_FRA{1,g}(idx_area,1) ;
        variable_area{1,g,3,area}(1:size(GNG_analysis.FR_coeff_var_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.FR_coeff_var_FRA {1,g}(idx_area,1);
        variable_area{1,g,4,area}(1:size(GNG_analysis.peak_latency_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.peak_latency_FRA{1,g}(idx_area,1) ;
        variable_area{1,g,5,area}(1:size(GNG_analysis.FWHM_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.FWHM_FRA{1,g}(idx_area,1);
        variable_area{1,g,6,area}(1:size(GNG_analysis.min_latency_FRA{1,g}(idx_area,1),1),1) = GNG_analysis.min_latency_FRA{1,g} (idx_area,1);
        variable_area{1,g,7,area}(1:size(GNG_analysis.Fraction_responsive_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.Fraction_responsive_FRA{1,g}(idx_area,1) ;
        variable_area{1,g,8,area}(1:size(GNG_analysis.lifetime_spar_FRA{1,g}(idx_area,1),1),1)  = GNG_analysis.lifetime_spar_FRA{1,g} (idx_area,1);

            for v = 1:size(variable_area,3)
                mean_var_area(1,g,v,area) = nanmean(variable_area{1,g,v,area}) ;
                size_n_area = size(variable_area{1,g,v,area}(~isnan(variable_area{1,g,v,area})),1) ;
                stderr_var_area(1,g,v,area) = nanstd(variable_area{1,g,v,area} ./sqrt(size_n_area)) ;
            end
        end   
        end
end




% compare them statistically in experts
for v = 1:size(variable_all,3)
     size_n_adol_all= size(variable_area{1,1,v}(~isnan(variable_area{1,1,v})),1) ;
        size_n_adult_all = size(variable_area{1,2,v}(~isnan(variable_area{1,2,v})),1) ;
        if ~isempty(size_n_adol_all)
            if ~isempty(size_n_adult_all)
                % compare the distribution 
                [p_all(v), h, stats] = ranksum(variable_all{1,1,v}, variable_all{1,2,v}, 'tail', 'both');

                effect = meanEffectSize(variable_all{1,1,v},variable_all{1,2,v},Paired=false,Effect="robustcohen",Alpha=0.05) ;

                CohenD_all(v) = effect.Effect ; 
                CI_all(:,v) = effect.ConfidenceIntervals ;

              
                % assign according to significancxe
                if  p_all(v) < 0.0005
                    p_txt{v} = '***';
                elseif  p_all(v) < 0.005
                    p_txt{v} = '**';
                elseif  p_all(v) < 0.05
                    p_txt{v} = '*';
                elseif p_all(v)> 0.05
                    p_txt{v} = 'n.s.';
                end
                
            end
        end


        % same per area
        for area = 1:length(areas)
            size_n_adol_area = size(variable_area{1,1,v,area}(~isnan(variable_area{1,1,v,area})),1) ;
            size_n_adult_area = size(variable_area{1,2,v,area}(~isnan(variable_area{1,2,v,area})),1) ;
            if size_n_adol_area > 0
                if size_n_adult_area > 0
                    [p_area(v,area), h, stats] = ranksum(variable_area{1,1,v,area}, variable_area{1,2,v,area}, 'tail', 'both');

                    effect = meanEffectSize(variable_area{1,1,v,area},variable_area{1,2,v,area},Paired=false,Effect="robustcohen",Alpha=0.05) ;

                    CohenD_area(v,area) = effect.Effect ;
                    CI_area(:,v,area) = effect.ConfidenceIntervals ;

                    if  p_area(v,area) < 0.0005
                        p_txt_area{v,area} = '***';
                    elseif p_area(v,area) < 0.005
                        p_txt_area{v,area} = '**';
                    elseif  p_area(v,area) < 0.05
                        p_txt_area{v,area} = '*';
                    elseif p_area(v,area)> 0.05
                        p_txt_area{v,area} = 'n.s.';
                    end
                elseif size_n_adult_area == 0
                     CohenD_area(v,area) = nan ; 
                    CI_area(:,v,area) = nan ;
                    p_area(v,area) = 1 ; 
                end
            elseif size_n_adol_area == 0
                 CohenD_area(v,area) = nan ; 
                CI_area(:,v,area) = nan  ;
                 p_area(v,area) = 1 ; 
 
            end

        end
end

clc
toc

%%
%% TABLE1 : table of average neuronal firing properties for all modulated neurons 
clc

% squeeze to 2D for  the overall mean and stderr
 for v = 1:size(variable_all,3)
        mean_var_2d_all(v,:) = squeeze (mean_var_all(1,:,v))' ;
        std_var_2d_all(v,:) = squeeze (stderr_var_all (1,:,v))' ;
    end
    
    mean_var_2d_all = round (mean_var_2d_all,4) ;
    std_var_2d_all = round (std_var_2d_all,4) ;
    CohenD_all = round(CohenD_all,4) ;
    CI_all = round(CI_all,4) ;
    
    GNG_ephys_table_all = table(categorical({...
        'spontaneous FR';...
        'evoked FR'; ...
        'FR coeff. var';...
        'latency to peak';...
        'FWHM';...
        'min. latency';...
        '% trials resp.';...
        'lifetime sparse.'}),...
        mean_var_2d_all(:,1),mean_var_2d_all(:,2),...
        std_var_2d_all(:,1),std_var_2d_all(:,2), CohenD_all',CI_all(1,:)',CI_all(2,:)',p_all'...
        ,'VariableNames',{'neuronal property',...
        'adol. mean',...
        'adult mean',...
        'adol. std',...
        'adult std',...
        'robust Cohen D',...
        'lower CI',...
        'upper CI',...
        'p-value'}) ;
   
    writetable(GNG_ephys_table_all,['GNG_ephys_table_exc_all_FRA.xlsx'],'Sheet',1) %; 

%% AREA Table 

for  area = 1:length(areas)
    if sum_mod_area_recs(area,1) > 5 & sum_mod_area_recs(area,2) > 5
        % squeeze to 2D for  the overall mean and stderr
        for v = 1:size(variable_all,3)
            mean_var_2d_area(v,:) = squeeze (mean_var_area(1,:,v,area))' ;
            std_var_2d_area(v,:) = squeeze (stderr_var_area (1,:,v,area))' ;
            CI_area_2d(v,:) = squeeze(CI_area(:,v,area)) ;
        end

        mean_var_2d_area = round (mean_var_2d_area,4) ;
        std_var_2d_area = round (std_var_2d_area,4) ;
        CI_area_2d = round (CI_area_2d,4) ;
        CohenD_area = round (CohenD_area,4) ;

        GNG_ephys_table_area = table(categorical({...
            'spontaneous FR';...
            'evoked FR'; ...
            'FR coeff. var';...
            'latency to peak';...
            'FWHM';...
            'min. latency';...
            '% trials resp.';...
            'lifetime sparse.'}),...
            mean_var_2d_area(:,1),mean_var_2d_area(:,2),...
            std_var_2d_area(:,1),std_var_2d_area(:,2),...
            CohenD_area(:,area), CI_area_2d(:,1), CI_area_2d(:,2), p_area(:,area)...
            ,'VariableNames',{'neuronal property',...
            'adol. mean',...
            'adult mean',...
            'adol. std',...
            'adult std',...
            'robust Cohen D',...
            'lower CI',...
            'upper CI',...
            'p-value'}) ;
        writetable(GNG_ephys_table_area,['GNG_ephys_table_area_FRA.xlsx'],'Sheet',area) %;
    end
end

%% table of area comparison 
for g = 1:4
       idx_area = sum_mod_area_recs(:,g) > 5
   

    end

for g = 1:numel(GNG_rec_all_cell)-2 % experts only

    idx_area = find(sum_mod_area_recs(:,g) > 5) ;

    clear p_all_area
    if length(idx_area) == 4

    for v = 1:size(variable_all,3)

        [~,~,stats] =  kruskalwallis([  squeeze(variable_area{1,g,v,1}),....
            squeeze(variable_area{1,g,v,2}),squeeze(variable_area{1,g,v,3}), squeeze(variable_area{1,g,v,4})],[],'off') ;
        stat_multi_all{g,v} = multcompare(stats,"Display","off");
        Ad_mean(:,v) =  squeeze(mean_var_area(1,g,v,1)) ;
        Ad_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,1)) ;
        Ap_mean(:,v) =  squeeze(mean_var_area(1,g,v,2)) ;
        Ap_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,2)) ;
        Av_mean(:,v) =  squeeze(mean_var_area(1,g,v,3)) ;
        Av_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,3)) ;
        TEa_mean(:,v) =  squeeze(mean_var_area(1,g,v,4)) ;
        TEa_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,4)) ;
        p_all_area(:,v) = stat_multi_all{g,v}(:,6) ;

    end

    GNG_ephys_table_group = table(categorical({...
        'spontaneous FR';...
        'evoked FR'; ...
        'FR coeff. var';...
        'latency to peak';...
        'FWHM';...
        'min. latency';...
        '% trials resp.';...
        'lifetime sparse.'}),...
        round(Ad_mean',4),  round(Ad_stderr',4), round(Ap_mean',4),  round(Ap_stderr',4),...
        round(Av_mean',4),  round(Av_stderr',4), round(TEa_mean',4),  round(TEa_stderr',4),...
        round(p_all_area(1,:)',4),  round(p_all_area(2,:)',4), round(p_all_area(3,:)',4), round(p_all_area(4,:)',4), round(p_all_area(5,:)',4), round(p_all_area(6,:)',4),...
        'VariableNames',{'neuronal property',...
        'Ad Mean',...
        'Ad STE',...
        'Ap Mean',...
        'Ap STE',...
        'Av Mean',...
        'Av STE',...
        'TEa Mean',...
        'TEa STE',...
        'AUDd - AUDp',...
        'AUDd - AUDv',...
        'AUDd - TEa',...
        'AUDp - AUDv',...
        'AUDp - TEa',...
        'AUDv - TEa',...
        }) ;
writetable(GNG_ephys_table_group,['GNG_ephys_table_group.xlsx'],'Sheet',g)
    



    elseif length(idx_area) == 3

    for v = 1:size(variable_all,3)

        [~,~,stats] =  kruskalwallis([  squeeze(variable_area{1,g,v,1}),....
            squeeze(variable_area{1,g,v,2}),squeeze(variable_area{1,g,v,3})],[],'off') ;
        stat_multi_all{g,v} = multcompare(stats,"Display","off");
        Ad_mean(:,v) =  squeeze(mean_var_area(1,g,v,1)) ;
        Ad_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,1)) ;
        Ap_mean(:,v) =  squeeze(mean_var_area(1,g,v,2)) ;
        Ap_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,2)) ;
        Av_mean(:,v) =  squeeze(mean_var_area(1,g,v,3)) ;
        Av_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,3)) ;
        p_all_area(:,v) = stat_multi_all{g,v}(:,6) ;

    end

    GNG_ephys_table_group = table(categorical({...
        'spontaneous FR';...
        'evoked FR'; ...
        'FR coeff. var';...
        'latency to peak';...
        'FWHM';...
        'min. latency';...
        '% trials resp.';...
        'lifetime sparse.'}),...
        round(Ad_mean',4),  round(Ad_stderr',4), round(Ap_mean',4),  round(Ap_stderr',4),...
        round(Av_mean',4),  round(Av_stderr',4),...
        round(p_all_area(1,:)',4),  round(p_all_area(2,:)',4), round(p_all_area(3,:)',4), round(p_all_area(4,:)',4), round(p_all_area(5,:)',4), round(p_all_area(6,:)',4),...
        'VariableNames',{'neuronal property',...
        'Ad Mean',...
        'Ad STE',...
        'Ap Mean',...
        'Ap STE',...
        'Av Mean',...
        'Av STE',...
        'AUDd - AUDp',...
        'AUDd - AUDv',...
        'AUDp - AUDv',...
 
        }) ;
writetable(GNG_ephys_table_group,['GNG_ephys_table_group.xlsx'],'Sheet',g)
    




    elseif length(idx_area) == 2



    for v = 1:size(variable_all,3)

        [~,~,stats] =  kruskalwallis([  squeeze(variable_area{1,g,v,2}),....
            squeeze(variable_area{1,g,v,3})],[],'off') ;
        stat_multi_all{g,v} = multcompare(stats,"Display","off");
       
        Ap_mean(:,v) =  squeeze(mean_var_area(1,g,v,2)) ;
        Ap_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,2)) ;
        Av_mean(:,v) =  squeeze(mean_var_area(1,g,v,3)) ;
        Av_stderr(:,v) =  squeeze(stderr_var_area(1,g,v,3)) ;
        p_all_area(:,v) = stat_multi_all{g,v}(:,6) ;

    end

    GNG_ephys_table_group = table(categorical({...
        'spontaneous FR';...
        'evoked FR'; ...
        'FR coeff. var';...
        'latency to peak';...
        'FWHM';...
        'min. latency';...
        '% trials resp.';...
        'lifetime sparse.'}),...
        round(Ap_mean',4),  round(Ap_stderr',4),...
        round(Av_mean',4),  round(Av_stderr',4),...
        round(p_all_area(1,:)',4),...
        'VariableNames',{'neuronal property',...
        'Ap Mean',...
        'Ap STE',...
        'Av Mean',...
        'Av STE',...       
        'AUDp - AUDv',...
        }) ;
writetable(GNG_ephys_table_group,['GNG_ephys_table_group.xlsx'],'Sheet',g)
    
    end 
end




%% FIGURES of all areas per group
clc
close all
v_label = ({...
    'spontaneous FR';...
    'evoked FR'; ...
    'FR coeff. var';...
    'latency to peak';...
    'FWHM';...
    'min. latency';...
    '% trials resp.';...
    'lifetime sparse.';})
Yscale_all = ({'linear';'linear';'linear';'linear';'linear';'linear';'linear';'linear';'linear'})

for g = 1:numel(GNG_rec_all_cell)-2 % experts only
    figure('units','normalized','outerposition',[0 0 1 1])
    for v = 1:size(variable_all,3)
        subplot(2,4,v)
        if sum_mod_area_recs(1,g) > 10 
            vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
            vec_length(:,1) = 1;
            boxchart(vec_length,squeeze(variable_area{1,g,v,1}),'BoxFaceColor',Colors_area{1},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
            hold on 
            errorbar(1, squeeze(mean_var_area(1,g,v,1)), squeeze(stderr_var_area(1,g,v,1)),'Color',Colors_area{1})
            plot(1, squeeze(mean_var_area(1,g,v,1)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{1},'MarkerEdgeColor',Colors_area{1})
            hold on
         

        end

        if sum_mod_area_recs(2,g) > 10 
            vec_length = nan(length(squeeze(variable_area{1,g,v,2})),1) ;
            vec_length(:,1) = 2;
            boxchart(vec_length,squeeze(variable_area{1,g,v,2}),'BoxFaceColor',Colors_area{2},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
            vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
            errorbar(2, squeeze(mean_var_area(1,g,v,2)), squeeze(stderr_var_area(1,g,v,2)),'Color',Colors_area{2})
            plot(2, squeeze(mean_var_area(1,g,v,2)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{2},'MarkerEdgeColor',Colors_area{2})
            hold on
        end

        if sum_mod_area_recs(3,g) > 10 
            vec_length = nan(length(squeeze(variable_area{1,g,v,3})),1) ;
            vec_length(:,1) = 3;
            boxchart(vec_length,squeeze(variable_area{1,g,v,3}),'BoxFaceColor',Colors_area{3},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
            errorbar(3, squeeze(mean_var_area(1,g,v,3)), squeeze(stderr_var_area(1,g,v,3)),'Color',Colors_area{3})
            plot(3, squeeze(mean_var_area(1,g,v,3)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{3},'MarkerEdgeColor',Colors_area{3})
            hold on
        end

        if sum_mod_area_recs(4,g) > 10 
            vec_length = nan(length(squeeze(variable_area{1,g,v,4})),1) ;
            vec_length(:,1) = 4;
            boxchart(vec_length,squeeze(variable_area{1,g,v,4}),'BoxFaceColor',Colors_area{4},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
            errorbar(4, squeeze(mean_var_area(1,g,v,4)), squeeze(stderr_var_area(1,g,v,4)),'Color',Colors_area{4})
            plot(4, squeeze(mean_var_area(1,g,v,4)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{4},'MarkerEdgeColor',Colors_area{4})

        end

        xticks([1 2 3 4])
        xticklabels({'AUDd','AUDp','AUDv','TEa'})
        box off
        ylabel(v_label(v))
        set(subplot(2,4, v),'YScale',Yscale_all(v));



    idx_area = find(sum_mod_area_recs(:,g) > 5) ;

    if length(idx_area) == 4
        max_y = [] ;
        max_y = max(max([  squeeze(variable_area{1,g,v,1}),....
            squeeze(variable_area{1,g,v,2}),squeeze(variable_area{1,g,v,3}), squeeze(variable_area{1,g,v,4})])) ;
        ylim([0 max_y*2])

        one_line = max_y + (max_y/5) ;
        two_line = max_y + ((max_y/5)*2) ;
        three_line = max_y + ((max_y/5)*3) ;
        four_line = max_y + (max_y/5) ;
        five_line = max_y + ((max_y/5)*4) ;
        six_line = max_y + (max_y/5) ;

        one_text = max_y + ((max_y/5) *1.5) ;
        two_text = max_y + ((max_y/5)*2.5) ;
        three_text = max_y + ((max_y/5)*3.5) ;
        four_text = max_y + ((max_y/5) *1.5)  ;
        five_text = max_y + ((max_y/5)*4.5) ;
        six_text = max_y + ((max_y/5) *1.5)  ;

        [~,~,stats] =  kruskalwallis([  squeeze(variable_area{1,g,v,1}),....
            squeeze(variable_area{1,g,v,2}),squeeze(variable_area{1,g,v,3}), squeeze(variable_area{1,g,v,4})],[],'off') ;
        stat_multi = multcompare(stats,"Display","off");
        p =stat_multi(:,6)

        line_left = [1.05 1.05 1.05 2.05 2.05 3.05] ;
        line_right =  [1.95 2.95 3.95 2.95 3.95 3.95] ;

        line_height = [one_line two_line three_line four_line five_line six_line] ;
        text_mean = [1.5 2 2.5 2.5 3 3.5] ;

        text_height = [one_text two_text three_text four_text five_text six_text] ;  ;

        for s = 1:6
            if  p(s) < 0.0005
                txt = '***';
            elseif p(s)  < 0.005
                txt = '**';
            elseif   p(s)  < 0.05
                txt = '*';
            elseif p(s)  > 0.05
                txt = 'ns';
            end

            line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
            hold on
            text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',10)

        end

        end
    end
end
directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';
filename = 'SuppFig412_adolescence'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');

filename = 'SuppFig412_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');

%% tuning bandwidth in oct. for AUDp and AUDv
clc
toc
close all
for area = 2:3
    if area == 2
        area_2 = 4; 
    elseif area == 3
        area_2 = 1 ; 
    end 
    idx_area = find(cell2mat(area_all_cell{1,g}) == area) ;
 idx_area_2 = find(cell2mat(area_all_cell{1,g}) == area_2) ;
 idx_area = [idx_area; idx_area_2] ;
 for g = 1:4
        vari{1,g}  = nan (150,1) ;
        vari{1,g}(idx_area) =  oct_mod_freq_FRA{1,g}(idx_area,:)
    end



figure
violinplot([vari{1,3} vari{1,1} vari{1,4} vari{1,2}]...
    ,{'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
ylim([0 4])
yticks([0:0.5:3])
ylabel('tuning bandwidth (oct.)')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ vari{1,3} vari{1,4} vari{1,1} vari{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95]
line_right =  [1.95 3.95 3.05 3.95]
line_height = [3.1 3.1 3.3 3.5]
text_mean = [1.5 3.5 2 3 ]
text_height = [3.2 3.2 3.4 3.6] 

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end




end 
%% population responsiveness for AUDp and AUDv
toc
for area = 2:3
    if area == 2
        area_2 = 1; 
    elseif area == 3
        area_2 = 1 ; 
    end 
    idx_area = find(cell2mat(area_all_cell{1,g}) == area) ;
 idx_area_2 = find(cell2mat(area_all_cell{1,g}) == area_2) ;
 idx_area = [idx_area; idx_area_2] ;
    for g = 1:4
        vari{1,g}  = nan (150,1) ;
        vari{1,g}(idx_area) =  sum_mod_freq_FRA{1,g}(idx_area,:)
    end


figure
violinplot([vari{1,3} vari{1,1} vari{1,4} vari{1,2}],...
    {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
ylim([0 1.6])
yticks([0:0.25:1])
ylabel('population responsiveness')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ vari{1,3} vari{1,4} vari{1,1} vari{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6) ;

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [1.1 1.1 1.3 1.5] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [1.2 1.2 1.4 1.6] ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
end

%% Bf distance to go 
close all
toc
clear vari
for area = 2:3
      if area == 2
        area_2 = 1; 
    elseif area == 3
        area_2 = 4 ; 
    end 
    idx_area = find(cell2mat(area_all_cell{1,g}) == area) ;
 idx_area_2 = find(cell2mat(area_all_cell{1,g}) == area_2) ;
 idx_area = sort([idx_area; idx_area_2]) ; 

 for g = 1:2
      vari{g,1}  = nan (150,1) ;

     vari{g,1}(idx_area) =  dist_easy_go{g,1}(idx_area,:)-0.1 ;
 end
 for g = 3:4
      vari{g,1}  = nan (150,1) ;

     vari{g,1}(idx_area) =  dist_easy_go{g,1}(idx_area,:) ;
 end


   clc
figure
violinplot([vari{3,1} vari{1,1} vari{4,1} vari{2,1}],...
       {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )

xlim([0 5])
%ylim([0 3.5])
%yticks([0:1:3])
ylabel('BF to Go (oct.)')
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;

 [~,~,stats] = kruskalwallis([ vari{3,1} vari{4,1} vari{1,1} vari{2,1}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6);

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [2.8 2.8 3 3.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [2.9 2.9 3.1 3.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
end
%% mean dprime for AUDp and AUDv
close all
clear vari mean_dprime_neuron
toc
for area = [2 3]
    if area == 2
        area_2 = 4; 
    elseif area == 3
        area_2 = 4 ; 
    end 
    idx_area = find(cell2mat(area_all_cell{1,g}) == area) ;
 idx_area_2 = find(cell2mat(area_all_cell{1,g}) == area_2) ;
 idx_area = [idx_area; idx_area_2] ;

for g   = 1:numel(GNG_rec_all_cell) % run per group
    mean_dprime_neuron{1,g} = nan(150,1) ;

    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        d_prime{1,g}(d_prime{1,g} == 0) = nan ; 
            mean_dprime_neuron{1,g}(c)  = nanmean(nanmean(squeeze(d_prime{1,g}(learned_stim,learned_stim,c)))) ;

    end


        vari{1,g}  = nan (150,1) ;
        vari{1,g}(idx_area) =  mean_dprime_neuron{1,g}(idx_area)
    
end 


   clc
          figure
violinplot([vari{1,3}  vari{1,1} vari{1,4}  vari{1,2} ],...
    {'Novice','Expert','Novice','Expert'},"ViolinColor",{[.5 .5 .5; .5 .5 .5; .0 .0 .0; .0 .0 .0]} )
xlim([0 5])
ylim([0 4.5])
yticks([0:1:4])
ylabel("freq. discrimination(d')")


ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
 [~,~,stats] = kruskalwallis([ vari{1,3} vari{1,4} vari{1,1} vari{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [3.8 3.8 4 4.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [3.9 3.9 4.1 4.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
end
