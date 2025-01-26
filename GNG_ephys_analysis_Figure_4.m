clear 
clc
close all
clearvars -except GNG_rec_all_cell
%%
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select GNG_rec_all_cell ');
addpath(path)
load (file)
cd (path)

%% Parameters
% area and layer parameters
areas= (1:4); % number of ACx subregions 
cortex = ({'AUDd6a', 'AUDp6a', 'AUDv6a', 'TEa6a',...
    'AUDd5', 'AUDp5', 'AUDv5', 'TEa5'}) ; % all subregions of ACx included in the analysis
area_str ={'AUDd','AUDp','AUDv','TEa'};
areas= 1:length(area_str) ; % number of ACx subregions 
area_str_reverse = {'TEa','AUDv','AUDp','AUDd'};
table_cat_areas = categorical({'AUDd'; 'AUDp'; 'AUDv'; 'TEa';})  ;

%group and stimulus parameters
group_str = {'adol.','adult','adol.','adult'};
n_learned_stimuli = 4; 

% time parameters
baseline = [-0.2 -0.05]; % baseline according to stimulus onset
start_baseline = 1; % baseline onset ms
stop_baseline = dist(baseline(1), baseline(2))* 1000; % baseline  ms
startStim = dist(baseline(1), 0)* 1000 ; % stimulus onset
stopStim = dist(baseline(1), 0.1)* 1000 ; % stimulus offset + 20 ms
window = [-0.2; 0.6]; % time window according to stimulus onset
startRange = window(1); %beginning of baseline
stopRange = window(2); %last timepoint before the reward
offset = 50 ; %offset after tone offset in ms 

% minimal criteria 
min_spikes = 100; % minimum spikes 

% iteration parameters
binSize = 0.001; % bin = 1 ms
rasterScale = 1; % binsize of raster plot
smoothSize = 5; % binsize of PSTh smoothing

% lick parameters
startRange_resp = -0.2 ;
stopRange_resp = 0.6 ;
binBorders = startRange_resp:binSize:stopRange_resp ; %define the time points
numBins = length(binBorders)-1 ;
startResp = 0 ; 
stopResp = 1.8 ;
startLick = -0.2 ;
stopLick = 0.2 ;

% colors
Colors_area = {[.1 .3 .8], [.5 .4 .9], [0 .5 .6],[0.4940 0.1840 0.5560]};
colors_ado_adu = {[.5 .7 .2],[ .2 .4 .2],[.5 .7 .2],[ .2 .4 .2]} ; 
L = {'-', '-','-','-'};
L2 = {'--', '-','--','-'};

directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';

%% enter the lion's den

tic
% test if neurons are modulated by the tone
for g   = 1:numel(GNG_rec_all_cell) % run per group

    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
     GNG_analysis.clusterIDs_all {1,g} = nan (300,length(Recs), length(areas)) ;
    GNG_analysis.clusterDepths_all{1,g} = nan (300,length(Recs), length(areas)) ;
    GNG_analysis.clusterIDs_mod{1,g} = nan(300,length(Recs), length(areas)) ;
    GNG_analysis.clusterDepths_mod{1,g} = nan(300,length(Recs), length(areas)) ;
    
    for i = 1:length(Recs) % run per recording

        spikes = GNG_rec_all_cell{1,g}(i).Units; %extract spikes
        spikeTimes = spikes.st; %extract spike times
        clusters = spikes.clu; %extract cluster ids per spike
        clusterIDs = unique(GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu); % extract the unique clusterIDs  
        
        eventTimes_all = GNG_rec_all_cell{1,g}(i).eventTimes.all;        

        %analyze the lick events
        eventTimes_licks = GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_licks ;
        eventTimes_reward = GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_reward  ;

        
        lick_index = zeros(1, length(eventTimes_licks)) ;
        
        for  Lick = 1:length (eventTimes_licks)
            L_smaller = eventTimes_licks(1,Lick) < min(eventTimes_all);
            L_within =  eventTimes_licks(1,Lick) > (eventTimes_all + stopResp)...
                & eventTimes_licks(1,Lick) < (eventTimes_all + startResp);
            L_larger = eventTimes_licks(1,Lick) > max(eventTimes_all);
            
            if ~isempty(L_smaller(L_smaller== 1))
                lick_index (1,Lick) = 1 ;
            elseif isempty(L_smaller(L_smaller== 1))
                lick_index (1,Lick) = 0 ;
            end
            
            if ~isempty(L_within(L_within== 1))
                lick_index (1,Lick) = 1 ;
            elseif isempty(L_within(L_within== 1))
                lick_index (1,Lick) = 0 ;
            end
            
            if ~isempty(L_larger(L_larger== 1))
                lick_index (1,Lick) = 1 ;
            elseif isempty(L_larger(L_larger== 1))
                lick_index (1,L_larger) = 0 ;
            end
        end
        
        eventTimes_licks = eventTimes_licks(lick_index == 1);
        if ~isempty(eventTimes_licks) 

        diff_lick_times = diff(eventTimes_licks) ;
       baseline_lick_idx =  find (diff_lick_times > 0.2) ;
       eventTimes_licks_baseline = eventTimes_licks(baseline_lick_idx) ;
        end
        
        
        eventTimes_licks_overview{g,i} = eventTimes_licks ;

        % run this for loop to determine significantly modulated units
        % across all events.
        for area = 1:length(areas)

            % choose only clusters in the relveant subregions and layers of ACx
                all_clusters = find(strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area)))...
                    |strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area+4)))) ;    
           clusterIDs = GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu(all_clusters) ; 
           clusterIDs =  clusterIDs (clusterIDs >0) ;


           % choose the depth of neruons 
           clusterDepths = GNG_rec_all_cell{1,g}(i).Units.Unit_table.spikeDepths(all_clusters);  %get spikedepth for a super sick graph later
           clusterDepths(clusterDepths>3850) = nan ; 
           clusterDepths = round((clusterDepths / 1000),2) ;            
           GNG_analysis.clusterIDs_all{1,g}(1:length(clusterIDs),i,area) = clusterIDs ;
           GNG_analysis.clusterDepths_all{1,g}(1:length(clusterDepths),i,area) = clusterDepths ; 

            for c = 1:length(clusterIDs)  %run per cluster = per neuron
                
                st = spikeTimes(clusters == clusterIDs(c));

                %sort out spike times outside of the event time range
                spikeTime = st(st>min(eventTimes_all+startRange) & st<max(eventTimes_all+stopRange));
                
                if length(spikeTime) > min_spikes
                    %extract the basics 
                    [~, psth,binArray,binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, eventTimes_all, spikeTime);
                    
                    % extract the smoothed basics
                    [~,psth_Smoothed, frArray_Smoothed,spikeCounts] ...
                        = GNG_smoothed_PSTH (startRange, stopRange, binCenters,binArray, binSize, smoothSize);

                    % extract the baseline of each spike 
                    baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                    onset_psth = psth(startStim + 1 :stopStim + offset) ; % onset + 50 ms after tone offset to account for delayed firing
 

                    %excitated neurons 
                    % compare the psth of tone and baseline 
                    [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                    
                    if  p < 0.05 % sort out non-excited units
                        % assign all relevant variables 
                        GNG_analysis.clusterIDs_mod{1,g}(c,i,area) = clusterIDs(c);
                        GNG_analysis.clusterDepths_mod{1,g}(c,i,area) = clusterDepths(c);
                        GNG_analysis.binArray_mod{1,g}{c,i,area} = binArray ;
                        GNG_analysis.binCenters{1,g}{c,i,area} = binCenters ; 
                        GNG_analysis.psth_Smoothed_mod{1,g}{c,i,area} = psth_Smoothed ;
                        GNG_analysis.frArray_Smoothed_mod{1,g}{c,i,area} = frArray_Smoothed ;
                        GNG_analysis.area_ID{1,g}(c,i,area) = areas(area) ;
                        GNG_analysis.spikeCounts{1,g}{c,i,area} = spikeCounts;
                    end

                    %look at both excited and inhibited neurons or population decoding  
                    [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','both');
                    if  p < 0.05 % sort out non-excited units
                        GNG_analysis.spiketimes_mod{1,g}{c,i,area} = st ;
                    end 

                end 
                
                % obtain units modulated by licks
                spikeTimeL = st(st>min(eventTimes_licks+startLick) & st<max(eventTimes_licks+stopLick));

                if ~isempty(eventTimes_licks_baseline) & length(eventTimes_licks_baseline) > 1

                    if length(spikeTimeL) > min_spikes

                        [~, psth_lick, binArray_lick, ~] = GNG_binning_PSTH (startLick,stopRange, binSize, eventTimes_licks, spikeTimeL);
                        [~, psth_lick_baseline,binArray_lick_baseline, ~] = GNG_binning_PSTH (startRange,stopRange, binSize, eventTimes_licks_baseline, st);

                        % extract the smoothed basics
                        [~,psth_Smoothed_lick, frArray_Smoothed_lick,spikeCounts] ...
                            = GNG_smoothed_PSTH (startRange, stopRange, binCenters,binArray_lick_baseline, binSize, smoothSize);

                        lick_baseline = psth_lick_baseline(start_baseline:stop_baseline);
                        lick_onset = psth_lick_baseline(startStim +1:stopStim + offset);


                        [p,~,~] = ranksum(lick_onset , lick_baseline,0.05,'tail','both');

                        if  p < 0.05
                            GNG_analysis.clusterIDs_mod{2,g}(c,i,area) = clusterIDs(c);
                            GNG_analysis.binArray_mod{2,g}{c,i,area} = binArray_lick ;
                            GNG_analysis.binArray_mod_base{2,g}{c,i,area} = binArray_lick_baseline ;
                            GNG_analysis.psth_Smoothed_mod{2,g}{c,i,area} = psth_Smoothed_lick ;
                            GNG_analysis.frArray_Smoothed_mod{2,g}{c,i,area} = frArray_Smoothed_lick ;
                        end
                    end
                end
            end  
        end
    end
end
toc

%% Table Supplemental 1-1 Number of modulated neruons per area
clc
% sum per group per area and per recording 

for g   = 1:numel(GNG_rec_all_cell) % run per group
    
   % get rec number
       Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    
    for i = 1:length(Recs) % run per recording
        
        % sum per area for ALL neurons 
        for area = 1:length(areas)
            idx_all_area = ~isnan(GNG_analysis.clusterIDs_all{1,g}(:,i,area))  ;
            sum_all_area(1,g,i,area) = size(idx_all_area(idx_all_area == 1),1) ;
        end
        % exclude nans
        idx_all = ~isnan(GNG_analysis.clusterIDs_all{1,g}(:,i,:)) ;
        sum_all(g,i) = size(idx_all(idx_all == 1),1) ;
        
        % sum per area for MODULATED neurons 
        for area = 1:length(areas)
            idx = find(GNG_analysis.clusterIDs_mod{1,g}(:,i,area)== 0) ;
            GNG_analysis.clusterIDs_mod{1,g}(idx,i,area) = nan ;
            idx_mod_area = ~isnan(GNG_analysis.clusterIDs_mod{1,g}(:,i,area))  ;
            sum_mod_area(1,g,i,area) = size(idx_mod_area(idx_mod_area == 1),1) ;
            idx = [] ;
        end       
        % exclude nans
        idx_mod = ~isnan(GNG_analysis.clusterIDs_mod{1,g}(:,i,:))  ;
        sum_mod(1,g,i) = size(idx_mod(idx_mod == 1),1) ;
    end
end

% sum per area per group    
for g   = 1:numel(GNG_rec_all_cell) % run per group
        format short 
        sum_mod_recs(g,:) = sum(squeeze(sum_mod(1,g,:))) ; % all modulated neurons in ACx
        sum_all_g(g,:) = sum(sum_all(g,:)) ; % all neurons in ACx
        frac_mod_all(g,:) =   round(sum_mod_recs(g,:)/ sum_all_g(g,:),2) ; % fraction of modulated units      

        sum_area(g,:,:) = squeeze(sum_all_area(1,g,:,:)) ; % per area all neurons
        sum_mod_area_all(g,:,:) =  squeeze(sum_mod_area(1,g,:,:)) ; % per area only modulated neurons
        
        for area = 1:length(areas)
            sum_area_recs(area,g) = sum(squeeze(sum_area(g,:,area))) ; % all units per area
            sum_mod_area_recs(area,g) = sum(squeeze(sum_mod_area_all(g,:,area))) ; % all modulated units per area
            frac_mod_area(area,g) =   round(sum_mod_area_recs(area,g)/ sum_area_recs(area,g),2) ; % fraction of modulated units per area       
        end
end

% create table for expert and naive recordings 
for g = 1:2:numel(GNG_rec_all_cell)-1
   % cd (path)
GNG_SU_area_table = table(table_cat_areas,sum_area_recs(:,g), sum_mod_area_recs(:,g), frac_mod_area(:,g),...
    sum_area_recs(:,g+1), sum_mod_area_recs(:,g+1), frac_mod_area(:,g+1),'VariableNames',...
    {'areas','neurons adol.','exc. neurons adol.','% adol.',...
        'neurons adults','exc. neurons adults','% adults'}) ;
writetable(GNG_SU_area_table,'GNG_ecx_area_table.xlsx','Sheet',g) ;
end 


% compare if the number of neurons is statistically different 
g1 = [1 1 2 3] ;
g2 = [2 3 4 4] ;
for g = 1:numel(GNG_rec_all_cell)
    for area = 1:length(areas)
        N1 = sum(sum_area_recs(area,g1(g))) ;
        n1 = sum(sum_mod_area_recs(area,g1(g))) ;
        N2 = sum(sum_area_recs(area,g2(g))) ;
        n2 = sum(sum_mod_area_recs(area,g2(g))) ;

        [tbl,chi2stat_distr(area,g),pval_distr(area,g)] = crosstab(...
            [ones(1, N1), 2*ones(1, N2)], ...
            [1*ones(1, n1), 2*ones(1, N1-n1), ...
            1*ones(1, n2), 2*ones(1, N2-n2)]) ;
    end
end
%% Figure 5 C & Supplementary Figure 4.1 Expert Single Units distrubtion per probe depth 

close all
clc
figure
for g   = 1:numel(GNG_rec_all_cell)-2 % only for experts 
    % retrieve all cluster depths
    cdepth_all = reshape(GNG_analysis.clusterDepths_mod{1,g}(:,:,:),[],1) ;
     idx_cdepth_all = find(~isnan(cdepth_all)) ;
     cdepth_all = cdepth_all(idx_cdepth_all) ;
     
     % plot them in a histogram 
     h =  histfit(cdepth_all,50) ;
     hold on 
     h(1).FaceColor = colors_ado_adu{g};
     h(2).Color = colors_ado_adu{g};
     h(1).EdgeColor = 'none' ;
     h(1).FaceAlpha = 0.4 ;
     h(2).LineStyle = L{g} ;
     h(2).LineWidth = 5 ;
       
     %makepretty;
   xlabel('spike depth (mm)')
    xlim([0 5])
    xticks([0:1:5])
    ylabel('n neurons')
     ylabel('n neurons')
     box off
     ax = gca;
     ax.XAxis.FontSize = 20;
     ax.YAxis.FontSize = 20;
     movegui('east');
    
     hold on
end 
  filename = 'Fig5c'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%%  supp Fig 51qb depth by area
for g   = 1:numel(GNG_rec_all_cell)-2 % only for experts 

    for area = 1:length(areas)
        cdepth_area = [] ;

        cdepth_area = reshape(GNG_analysis.clusterDepths_mod{1,g}(:,:,area),[],1) ;
        idx_cdepth = find(~isnan(cdepth_area)) ;
        cdepth_area = cdepth_area(idx_cdepth) ;
        
       f =  figure(g+1) ;
       h =  histfit(cdepth_area,50) ;
       h(1).FaceColor = Colors_area{area};
       h(1).FaceAlpha = 0.5 ;
       h(2).Color = Colors_area{area};
       h(1).EdgeColor = 'none' ;
       h(1).EdgeColor = 'none' ;
       ylim([0 20 ])
       hold on 
    end

    makepretty;
    xlabel('spike depth (mm)')
    xlim([0 5])
    xticks([0:1:5])
    ylabel('n neurons')
    title(group_str{g})
    box off
    
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    
end

  filename = 'suppFig51b_adolescent'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
filename = 'suppFig51b_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(3);
saveas(recent_figure, filepath, 'svg');
%% Supplementary Figure 5.1 C Naive Single Units distrubtion per probe depth 
close all
clc

f1 =  figure(1) ;
for g   = 3:numel(GNG_rec_all_cell) % only for naive 
    % retrieve all cluster depths
    cdepth_all = reshape(GNG_analysis.clusterDepths_mod{1,g}(:,:,:),[],1) ;
     idx_cdepth_all = find(~isnan(cdepth_all)) ;
     cdepth_all = cdepth_all(idx_cdepth_all) ;
     
     figure(1)
     % plot them in a histogram 
     h =  histfit(cdepth_all,50) ;
     hold on 
     h(1).FaceColor = colors_ado_adu{g};
     h(2).Color = colors_ado_adu{g};
     h(1).EdgeColor = 'none' ;
     h(1).FaceAlpha = 0.4
     h(2).LineStyle = L{g} ;
     h(2).LineWidth = 5 ;
       
     %makepretty;
     xlabel('spike depth (\mum)')
     xlim([0 5])
     ylim
     xticks([0:1:5])
     xticklabels([0:1000:5000])
     ylabel('n neurons')
     box off
     ax = gca;
     ax.XAxis.FontSize = 20;
     ax.YAxis.FontSize = 20;
     movegui('east');
end

  filename = 'Fig7b'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');

%% Fig 5d psth by depth of example rec 
close all
clc
g = 2
i = 5
depthBinSize = 50; % in units of the channel coordinates, in this case ??m
timeBinSize = 0.01; % seconds
bslWin = [-0.2 -0.1]; % window in which to compute "baseline" rates for normalization
probe_length = 3850; % ms 
lower_limit = 1000 ; % most superficial probe limit to plot 
depthBins = lower_limit:depthBinSize:probe_length; % all bins 
nD = length(depthBins)-1;% length of depthbins 
window = [-0.2; 0.6]; % time window according to stimulus onset

% extract spike depths 
        spikeTimes = GNG_rec_all_cell{1,g}(i).Units.st ;
        spikeDepths = GNG_rec_all_cell{1,g}(i).Units.spikeDepths ;
        spikeDepths = probe_length - spikeDepths ;
        eventTimes = GNG_rec_all_cell{1,g}(i).eventTimes.all;
        
        % run per depth bin 
        for d = 1:nD
            theseSp = spikeDepths>depthBins(d) & spikeDepths<=depthBins(d+1);
            %extract psth per depth 
            [psth, timeBins, ~, ~, ~, ~] = psthAndBA(spikeTimes(theseSp), eventTimes, window, timeBinSize);
            
            if d==1
                allP = zeros(nD, length(psth));
            end
            allP(d,:) = psth;
        end

        figure;
        hold on 
        nD = length(depthBins)-1;
        ax = gca;
        imagesc(timeBins, depthBins(1:end-1), mat2gray(allP));
        hold on
        set(gca, 'YDir', 'reverse');
        hold on;
        xline(0,'w--')
        xline(0.1,'w--')
        xlabel(['time(sec)']);
        xlim([-0.2 0.6])
        ylim([1500 3800])
        ylabel('depth')
        box off
        h = colorbar
        h.Label.String = 'norm FR';
        makepretty
        box off
        ax = gca;
        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 10;
        movegui('east');

  filename = 'Fig5d_psth_bydepth'
filepath = fullfile(directory, filename);
recent_figure = gcf;
saveas(recent_figure, filepath, 'svg');
%% Figure 4 & Figure 6 population  PSTH all events  calculation 
% calculate the psth, fr array and bin array, as well as spiketimes per
% group for all neurons per AREA 
clc
for   g = 1:numel(GNG_rec_all_cell) % run per group

    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    
    psth_all_cell = GNG_analysis.psth_Smoothed_mod{1,g}(~cellfun('isempty', GNG_analysis.psth_Smoothed_mod{1,g})) ;
    psth_all_mat{1,g} = cell2mat(psth_all_cell) ;
    
    mean_pop_psth_all(g,:) = nanmean(psth_all_mat{1,g}) ;
    stderr_pop_psth_all(g,:) =  nanstd(psth_all_mat{1,g})./sqrt(size(psth_all_mat{1,g},1)) ;
    
    for area = 1:length(areas)
        for i = 1:length(Recs)
            for c = 1:size(GNG_analysis.psth_Smoothed_mod{1,g},1)
                psth_all_cell_areas{g,area}{c,i} = GNG_analysis.psth_Smoothed_mod{1,g}{c,i,area} ;
                
                spiketimes_all_cell_areas{g,area}{c,i} = GNG_analysis.spiketimes_mod{1,g}{c,i,area} ;
                spiketimes_all_cell_areas_i{g,i,area}{c} = GNG_analysis.spiketimes_mod{1,g}{c,i,area} ;
                
            end
            
            spiketimes_all_cell_areas_i{g,i,area} =  spiketimes_all_cell_areas_i{g,i,area}(~cellfun('isempty',  spiketimes_all_cell_areas_i{g,i,area})) ;
            spiketimes_all_cell_g{g}{i,area} =  spiketimes_all_cell_areas_i{g,i,area}(~cellfun('isempty',  spiketimes_all_cell_areas_i{g,i,area})) ;

        end
        psth_area_cell =  psth_all_cell_areas{g,area}(~cellfun('isempty',  psth_all_cell_areas{g,area})) ;
        psth_area_mat = cell2mat(psth_area_cell);
        
        mean_pop_psth_area{g}(area,:) = nanmean(psth_area_mat) ;
        stderr_pop_psth_area{g}(area,:) =  nanstd(psth_area_mat)./sqrt(size(psth_area_mat,1)) ;
        
        spiketimes_area_cell{g,area} = spiketimes_all_cell_areas{g,area}(~cellfun('isempty',  spiketimes_all_cell_areas{g,area})) ;

    end
end       
 
%% for Figure 5 Population decoding extract the raw spike times per recording 
%per recording all 
% calculate the psth, fr array and bin array, as well as spiketimes per group for all neurons
% group for all neurons per RECORDING 
clc
for   g = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    for i = 1:length(Recs)
        for area = 1:length(areas)
            for c = 1:size(GNG_analysis.psth_Smoothed_mod{1,g},1)
                spiketimes_all_cell_i{g,i}{c,area} = GNG_analysis.spiketimes_mod{1,g}{c,i,area} ;
            end
        end
        spiketimes_i_cell{g,i} =  spiketimes_all_cell_i{g,i}(~cellfun('isempty',  spiketimes_all_cell_i{g,i})) ;
        eventtimes_i_cell{g,i} = GNG_rec_all_cell{1, g}(i).eventTimes.all  ;
    end
end

clear st 
% extract spiketimes and eventtimes for population decoding in python
for   g = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    for i = 1:length(Recs)
        for c = 1:size(GNG_analysis.psth_Smoothed_mod{1,g},1)

            st_Ad{g,i}{c,1} = spiketimes_all_cell_i{g,i}{c,1} ;
            st_A1{g,i}{c,1} = spiketimes_all_cell_i{g,i}{c,2} ;
            st_Av{g,i}{c,1} = spiketimes_all_cell_i{g,i}{c,3} ;
            st_TEa{g,i}{c,1} = spiketimes_all_cell_i{g,i}{c,4} ;
        end

        st_Ad{g,i} =  st_Ad{g,i}(~cellfun('isempty',  st_Ad{g,i})) ;
        st_A1{g,i} =  st_A1{g,i}(~cellfun('isempty',  st_A1{g,i})) ;
        st_Av{g,i} =  st_Av{g,i}(~cellfun('isempty',  st_Av{g,i})) ;
        st_TEa{g,i} =  st_TEa{g,i}(~cellfun('isempty',  st_TEa{g,i})) ;

    end
end

st = spiketimes_i_cell ; 
et = eventtimes_i_cell;
%% calculate the Lick PSTH & correlation to neuronal PSTH
tic
clc

for g   = 1:numel(GNG_rec_all_cell) % run per group

    for rec = 1:numel(GNG_rec_all_cell{1,g})
        GNGtable_cell{g,rec} = GNG_rec_all_cell{1, g}(rec).Behavior;
    end

    for i = 1:numel(GNG_rec_all_cell{1,g})

        eventTimes_all = GNG_rec_all_cell{1,g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_licks ;

        if ~isempty(GNGtable_cell{g,i})

            if isa(GNGtable_cell{g,i}.stim_types, 'double')
                stim_ids = GNGtable_cell{g,i}.stim_ids;
                lick_times = GNGtable_cell{g,i}.lick_times;
                stim_times = GNGtable_cell{g,i}.stim_times;

            elseif isa(GNGtable_cell{g,i}.stim_types, 'cell')
                stim_ids  = GNGtable_cell{g,i}.stim_ids{1,1};
                lick_times = GNGtable_cell{g,i}.lick_times{1,1};
                stim_times = GNGtable_cell{g,i}.stim_times{1,1};
            end

        end
        %lick psth
        binArray = zeros(length(stim_times), numBins) ;

        for r = 1:length(stim_times)
            [n,binCenters] = histdiff(lick_times, stim_times(r), binBorders); %get the lick per binsize
            binArray(r,:) = n;
        end

        inclRange = binCenters > startRange_resp & binCenters<= stopRange_resp;
        lickCounts = sum(binArray(:,inclRange),2)./(stopRange_resp - startRange_resp); % sum the spikecounts

        gaussian_window = gausswin(round(smoothSize*6),3); %apply a gaussian window
        summed_window = gaussian_window./sum(gaussian_window); %sum all winows

        binArray_smoothed = conv2(summed_window,1,binArray', 'same')'./binSize;  %conv the summed window
        frArray_lick_all{g,i} = binArray_smoothed ;
        psth_lick_Smoothed(g,i,:) = mean(binArray_smoothed) ;

    
    end
end

%% Figure 4 c and supplementary Figure 6.1 D pop psth + lick psth
clc
close all
for g = 1:numel(GNG_rec_all_cell)
    npsth_Smoothed_all = [] ;
    nlick_Smoothed_all = [] ; 
     npsth_Smoothed_all = [] ; 
if g ==1 
    figure
elseif g == 3
    figure
end 

   for c = 1:size(psth_all_mat{1,g},1)
       psth_all_mat{1,g}(c,:) =  zscore(psth_all_mat{1,g}(c,1:end)) ;
       mean_pop_psth = mean (psth_all_mat{1,g}(c,1:startStim)) ;
       npsth_Smoothed_all(c,:) = abs(mean_pop_psth - psth_all_mat{1,g}(c,1:end)) ;
   end
    mean_psth_Smoothed_all = mean(npsth_Smoothed_all) ;
    stde_psth_Smoothed_all = std(npsth_Smoothed_all)./sqrt(size(psth_all_mat{1,g},1)) ;

        for i = 1:numel(GNG_rec_all_cell{1,g})
       psth_lick_Smoothed(g,i,:) =  zscore(psth_lick_Smoothed(g,i,1:end)) ;
       mean_pop_psth = mean (psth_lick_Smoothed(g,i,1:startStim)) ;
       nlick_Smoothed_all(i,:) = squeeze(abs(mean_pop_psth - psth_lick_Smoothed(g,i,1:end))) ;
     end
   
     mean_lick_Smoothed_all = mean (nlick_Smoothed_all ) ;
     stde_lick_Smoothed_all =  std(nlick_Smoothed_all)./sqrt(size(nlick_Smoothed_all,1)) ;
     
     mean_pop_psth = mean (mean_psth_Smoothed_all(1,1:startStim)) ;
     psth_Smoothed_all = abs(mean_pop_psth - mean_psth_Smoothed_all(1,1:end)) ;
     
     mean_pop_psth = mean (mean_lick_Smoothed_all(1,1:startStim)) ;
     lick_Smoothed_all = abs(mean_pop_psth - mean_lick_Smoothed_all(1,1:end)) ;

    plot([1,1:(dist(window(1), window(2))*1000)-1],...
        psth_Smoothed_all,'Color',colors_ado_adu{g},...
        'linestyle','-','linewidth',3);
    hold on
    patch([[1,1:(dist(window(1), window(2))*1000-1)]...
        flip([1,1:(dist(window(1), window(2))*1000-1)])],...
        [psth_Smoothed_all + ....
        stde_psth_Smoothed_all ...
        flip(psth_Smoothed_all -  stde_psth_Smoothed_all)]...
        , colors_ado_adu{g} ,...
        'facealpha' , 0.3, 'EdgeColor','none')
 
    hold on 
    
 plot([1,1:(dist(window(1), window(2))*1000)-1],...
        lick_Smoothed_all,'Color',colors_ado_adu{g},...
        'linestyle','--','linewidth',3);
    hold on
    patch([[1,1:(dist(window(1), window(2))*1000-1)]...
        flip([1,1:(dist(window(1), window(2))*1000-1)])],...
        [lick_Smoothed_all + ....
        stde_lick_Smoothed_all ...
        flip(lick_Smoothed_all -  stde_lick_Smoothed_all)]...
        , colors_ado_adu{g} ,...
        'facealpha' , 0.3, 'EdgeColor','none')
      
    xline(200,'--','Color','k');
    xline(300,'--','Color','k');
    xlabel('time(sec)')
    xlim([10 790])
    xticks([10 200 400 600 790])
    xticklabels([-200 0 200 400 600])
    ylim([-0.5 4])
    ylabel('normalized FR & LR')
    box off;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
lick_FR_corr{g} = corrcoef (lick_Smoothed_all,psth_Smoothed_all) ;
end





 filename = 'Fig5d'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
 filename = 'Fig7c'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');

%% Supplementary Figure 4.1 B plot population PSTh per auditory subregion 
close all 
for area = 1:length(areas)
    figure
    for g = 1:numel(GNG_rec_all_cell)-2 % only experts

        plot([1,1:(dist(window(1), window(2))*1000)-1],...
            mean_pop_psth_area{g}(area,:),'Color',Colors_area{area},...
            'linestyle',L2{g},'linewidth',2);
        hold on
        patch([[1,1:(dist(window(1), window(2))*1000-1)]...
            flip([1,1:(dist(window(1), window(2))*1000-1)])],...
            [mean_pop_psth_area{g}(area,:) + ....
            stderr_pop_psth_area{g}(area,:) ...
            flip(mean_pop_psth_area{g}(area,:) -  stderr_pop_psth_area{g}(area,:))]...
            , Colors_area{area} ,...
            'facealpha' , 0.2, 'EdgeColor','none')

        xline(200,'--','Color','k');
        xline(300,'--','Color','k');
        xlabel('time(sec)')
        xticks(0:200:800)
        xticklabels([-200:200:600])
        %ylim([-2 6])
        ylabel('normalized FR (Hz)')
        box off;
        hold on;
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('west');

    end
end

 filename = 'suppFig41c_AUDd'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');
 filename = 'suppFig41c_AUDp'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
 filename = 'suppFig41c_AUDv'
filepath = fullfile(directory, filename);
recent_figure = figure(3);
saveas(recent_figure, filepath, 'svg');
 filename = 'suppFig41c_TEa'
filepath = fullfile(directory, filename);
recent_figure = figure(4);
saveas(recent_figure, filepath, 'svg');
%% Table 1 - calculate the neuronal firing properties per best frequency
clc
% BF of modulated neurons in task
%allocate 
GNG_analysis.max_FR_disc = cell(1,2) ;
GNG_analysis.BF_idx_disc = cell(1,2) ;

tic
% extract BF and Max FR per cell per recording per area per group 
for g   = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;

    GNG_analysis.max_FR_disc{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.norm_max_FR_disc{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.BF_idx_disc{1,g} = nan (100,length(Recs), length(areas)) ;
    
    for i = 1:length(Recs) % run per recording
        for area = 1:length(areas)
            frarray_rec = {} ;
            frarray_rec = GNG_analysis.frArray_Smoothed_mod{1,g}(:,i,area) ;
            frarray_rec = frarray_rec(~cellfun('isempty', frarray_rec)) ;
            
            for c = 1:length(frarray_rec)  %run per cluster = per neuron
                psth = [] ; 
                FR_peaks= []; 
                max_stim = [] ; 
                
                for e = [5 6 11 12] % easy and hard go and no go 
                    % extract the relevant trials
                    idx_e =  GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx(e,:);
                    idx_e(isnan(idx_e)) = [] ;
                    idx_e = idx_e(1:end-1) ;
                    
                    psth(e,:) = mean(frarray_rec{c,1}(idx_e,:)) ;
                    
                    FR_peaks = findpeaks(psth(e,startStim +1 :stopStim + 50)) ; % if neuron has multiple peaks 
                    if ~isempty(FR_peaks)
                        max_FR =  max(FR_peaks) ; % the peak of peaks 
                    elseif isempty(FR_peaks)
                        max_FR = 0 ;
                    end
                    max_stim(e,:) = max_FR ;
                end
                
                [max_FR_disc, ~] = max(max_stim) ;
                
                BF_idx_disc_all =  find(max_stim == max_FR_disc) ; 
                
                if size(BF_idx_disc_all,1) > 1 % if the max FR of a neuron is through mutliple FR it fires non-preferential 
                    BF_idx_disc = 0 ;
                elseif size(BF_idx_disc_all,1) == 1
                    BF_idx_disc = BF_idx_disc_all ;
                end
                                
                GNG_analysis.max_FR_disc{1,g}(c,i,area) = max_FR_disc ;
                GNG_analysis.BF_idx_disc{1,g}(c,i,area) = BF_idx_disc ;
            end
        end
    end
end



% extract the basic neuronal features according to the best Freq of the neurons 
clc
%allocate 
GNG.analysis.baseline_FR = cell (1,2) ;
GNG_analysis.max_FR = cell(1,2) ; 
GNG_analysis.peak_latency = cell(1,2) ; 
GNG_analysis.Fraction_responsive = cell(1,2) ; 
GNG_analysis.FR_coeff_var = cell(1,2) ; 
GNG_analysis.min_latency = cell(1,2) ; 
GNG_analysis.lifetime_spar = cell(1,2) ;
GNG_analysis.FWHM = cell(1,2) ;


% neuronal properties per group per area per rec per cell 
for g   = 1:numel(GNG_rec_all_cell) % run per group
    
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;

    GNG_analysis.baseline_FR{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.max_FR{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.peak_latency{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.Fraction_responsive{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.FWHM{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.FR_coeff_var{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.min_latency{1,g} = nan (100,length(Recs), length(areas)) ;
    GNG_analysis.lifetime_spar{1,g} = nan (100,length(Recs), length(areas)) ;
    

    for i = 1:length(Recs) % run per recording
        
        for area = 1:length(areas)
           
            frarray_rec = {} ;
            frarray_rec = GNG_analysis.frArray_Smoothed_mod{1,g}(:,i,area) ;
            frarray_rec = frarray_rec(~cellfun('isempty', frarray_rec)) ;
            
            binArray_rec = {} ;
            binArray_rec = GNG_analysis.binArray_mod{1,g}(:,i,area) ;
            binArray_rec = binArray_rec(~cellfun('isempty', binArray_rec)) ;
            
            clustermod = GNG_analysis.clusterIDs_mod{1, g}(:,i,area)
            clustermod (clustermod == 0) = [] ; 


            for c = 1:length (frarray_rec) %run per cluster = per neuron
                
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
                jitter = [] ;
                min_latency =[] ;
                std_interval = [] ;
                neuronal_jitter = [];
                normalized_jitter = [] ;
                frarray_window_all  = [] ;
                fr_trial_all =[] ;
                lifetime_spar = [] ;
                norm_lifetime_spar = [] ;
               
                
                if GNG_analysis.BF_idx_disc{1,g}(c,i,area) == 0
                    BF_e = GNG_rec_all_cell{1,g}(i).eventTimes.tr_resp_idx(:,:) ;
                elseif GNG_analysis.BF_idx_disc{1,g}(c,i,area) ~= 0
                    BF_c = GNG_analysis.BF_idx_disc{1,g}(c,i,area) ;
                    BF_e = GNG_rec_all_cell{1,g}(i).eventTimes.tr_resp_idx(BF_c,:) ;
                end
                
                BF_e = BF_e(~isnan(BF_e)) ;
                BF_e = BF_e(:,1:end-1) ;
                
                psth = mean(frarray_rec{c,1}(BF_e,:)) ;
                baseline_FR = psth(start_baseline:stop_baseline) ;
                mean_baseline_FR =  mean(baseline_FR);
                GNG_analysis.baseline_FR{1,g}(c,i,area) = mean_baseline_FR ;
                
                psth_window = psth(startStim +1 :stopStim + 50) ;

                FR_peaks_stim_window = findpeaks(psth_window) ;
                FR_peaks_whole_window = findpeaks(psth(stop_baseline+50:end)) ;
                
                if ~isempty(FR_peaks_stim_window)
                    max_FR =  max(FR_peaks_stim_window) ;
                    GNG_analysis.max_FR{1,g}(c,i,area) = max_FR;
                    
                    % calculate the peak latency
                    max_FR_whole_window =  max(FR_peaks_whole_window) ;
                    index_peak_latency = find(psth(stop_baseline+50:end)  == max_FR_whole_window) ; % ms
                    peak_latency = max(index_peak_latency) ;
                    GNG_analysis.peak_latency{1,g}(c,i,area) = peak_latency ;

                   idx =  find(GNG_rec_all_cell{1, g}(i).Units.Unit_table.clu == clustermod(c))
                    GNG_rec_all_cell{1, g}(i).Units.Unit_table.peak_latency(idx) = peak_latency ;

                    %FWHM
                    half_max = max_FR / 2;
                    left_idx = find( psth(stop_baseline+50:end) > half_max, 1, 'first');
                    right_idx = find( psth(stop_baseline+50:end) > half_max, 1, 'last');
                    FWHM = (right_idx - left_idx);
                    GNG_analysis.FWHM{1,g}(c,i,area) = FWHM ;
                    GNG_rec_all_cell{1, g}(i).Units.Unit_table.FWHM(idx) = FWHM ;

                    
                    frarray_window  =  frarray_rec{c,1}(BF_e,startStim +1:stopStim + 50); % tone onset:tone offset + 50 ms
                    fr_trial = mean(frarray_window');
                    
                    % Fraction of responsive trials per neuron
                    nonactive_trials = find (fr_trial == 0) ;
                    n_trials = size(fr_trial,2) ;
                    Fraction_responsive = 1 - (length(nonactive_trials)/n_trials) ;
                    GNG_analysis.Fraction_responsive{1,g}(c,i,area) = Fraction_responsive ;
                    
                    % FR coefficient of variation
                    var_fr = var(fr_trial) ;
                    mean_fr = mean(fr_trial) ;
                    FR_coeff_var = sqrt(var_fr) / mean_fr;
                    GNG_analysis.FR_coeff_var{1,g}(c,i,area) = FR_coeff_var ;

                    %binarray_window = binArray_rec{c,1}(BF_e,startStim +1:stopStim + 50) ;                    
                    % extract the latency per trial (first
                    % spike and the spike iti 
                    
                    binarray_window = frarray_rec{c,1}(BF_e,startStim +1:stopStim + 50) ;
                    min_latency_trial = [] ;
                    spike_iti= [];
                 
                    for t = 1:size(binArray_rec{c,1}(BF_e,:),1)
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
                    GNG_analysis.min_latency{1,g}(c,i,area) = min_latency ;
                    GNG_rec_all_cell{1, g}(i).Units.Unit_table.min_latency(idx) = min_latency ;


                    % all learned tones 
                    fr_trial_stim = nan (13,1) ; 
                    for e = [5 6 11 12]
                        idx_e =  GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx(e,:);
                        idx_e(isnan(idx_e)) = [] ;
                        idx_e = idx_e(1:end-1) ;
                        frarray_window_stim = (frarray_rec{c,1}(idx_e,startStim+1: stopStim+50)) ;
                        fr_trial_stim(e,:) = nanmean(nanmean(frarray_window_stim'));
                    end 
                    
                    fr_trial_stim = fr_trial_stim(~isnan(fr_trial_stim)) ; 
                    % extract life time sparseness 
                    lifetime_spar = ( 1 - ( ( (sum( fr_trial_stim ) / n_learned_stimuli )^2) /...
                        ( (sum( fr_trial_stim.^2)) /n_learned_stimuli) ) )/( 1 - ( 1 / n_learned_stimuli ) )   ; 
                    % divide by the 4 different learned stimuli
                   % norm_lifetime_spar = abs(lifetime_spar)./100 ;
                    GNG_analysis.lifetime_spar{1,g}(c,i,area) = abs(lifetime_spar) ;
                    
                end              
            end
        end
    end
end










clear variable mean_var 
% calculate the population mean and std error of all neuronal features 
for g   = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
% extract all features 
        variable_all{1,g,1} = reshape (GNG_analysis.baseline_FR{1,g}(:,:,:),[],1) ;
        variable_all{1,g,2}  = reshape (GNG_analysis.max_FR{1,g}(:,:,:),[],1) ;
        variable_all{1,g,3}  = reshape (GNG_analysis.FR_coeff_var{1,g}(:,:,:),[],1) ;
        variable_all{1,g,4}  = reshape (GNG_analysis.peak_latency{1,g}(:,:,:),[],1) ;
        variable_all{1,g,5}  = reshape (GNG_analysis.FWHM{1,g}(:,:,:),[],1) ;
        variable_all{1,g,6} = reshape (GNG_analysis.min_latency{1,g}(:,:,:),[],1) ;
        variable_all{1,g,7}  = reshape (GNG_analysis.Fraction_responsive{1,g}(:,:,:),[],1) ;
        variable_all{1,g,8}  = reshape (GNG_analysis.lifetime_spar{1,g}(:,:,:),[],1) ;
        
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
            variable_area{1,g,1,area} = reshape (GNG_analysis.baseline_FR{1,g}(:,:,area),[],1) ;
            variable_area{1,g,2,area} = reshape (GNG_analysis.max_FR{1,g}(:,:,area),[],1) ;
            variable_area{1,g,3,area} = reshape (GNG_analysis.FR_coeff_var{1,g}(:,:,area),[],1) ;
            variable_area{1,g,4,area} = reshape (GNG_analysis.peak_latency{1,g}(:,:,area),[],1) ;
            variable_area{1,g,5,area} = reshape (GNG_analysis.FWHM{1,g}(:,:,area),[],1) ;
            variable_area{1,g,6,area} = reshape (GNG_analysis.min_latency{1,g}(:,:,area),[],1) ;
            variable_area{1,g,7,area} = reshape (GNG_analysis.Fraction_responsive{1,g}(:,:,area),[],1) ;
            variable_area{1,g,8,area} = reshape (GNG_analysis.lifetime_spar{1,g}(:,:,area),[],1) ;  
            
            for v = 1:size(variable_area,3)
                mean_var_area(1,g,v,area) = nanmean(variable_area{1,g,v,area}) ;
                size_n_area = size(variable_area{1,g,v,area}(~isnan(variable_area{1,g,v,area})),1) ;
                stderr_var_area(1,g,v,area) = nanstd(variable_area{1,g,v,area} ./sqrt(size_n_area)) ;
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
            if ~isempty(size_n_adol_area)
                if ~isempty(size_n_adult_area)
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
                    
                end
            end
        end
end
clc
toc

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
   
    writetable(GNG_ephys_table_all,['GNG_ephys_table_exc_all.xlsx'],'Sheet',1) %; 

%% AREA Table 

for  area = 1:length(areas)
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
writetable(GNG_ephys_table_area,['GNG_ephys_table_area.xlsx'],'Sheet',area) %;
end


%% table of area comparison 

for g = 1:numel(GNG_rec_all_cell)-2 % experts only
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
        vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
        vec_length(:,1) = 1; 
        boxchart(vec_length,squeeze(variable_area{1,g,v,1}),'BoxFaceColor',Colors_area{1},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
        hold on
        errorbar(1, squeeze(mean_var_area(1,g,v,1)), squeeze(stderr_var_area(1,g,v,1)),'Color',Colors_area{1})
        plot(1, squeeze(mean_var_area(1,g,v,1)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{1},'MarkerEdgeColor',Colors_area{1})
        vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
        vec_length(:,1) = 2; 
        boxchart(vec_length,squeeze(variable_area{1,g,v,2}),'BoxFaceColor',Colors_area{2},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
        vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
            errorbar(2, squeeze(mean_var_area(1,g,v,2)), squeeze(stderr_var_area(1,g,v,2)),'Color',Colors_area{2})
        plot(2, squeeze(mean_var_area(1,g,v,2)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{2},'MarkerEdgeColor',Colors_area{2})
        hold on 
        vec_length(:,1) = 3; 
        boxchart(vec_length,squeeze(variable_area{1,g,v,3}),'BoxFaceColor',Colors_area{3},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
                    errorbar(3, squeeze(mean_var_area(1,g,v,3)), squeeze(stderr_var_area(1,g,v,3)),'Color',Colors_area{3})
        plot(3, squeeze(mean_var_area(1,g,v,3)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{3},'MarkerEdgeColor',Colors_area{3})
        hold on
        vec_length = nan(length(squeeze(variable_area{1,g,v,1})),1) ;
        vec_length(:,1) = 4; 
        boxchart(vec_length,squeeze(variable_area{1,g,v,4}),'BoxFaceColor',Colors_area{4},'MarkerStyle','.','MarkerColor',[0.5 0.5 0.5])
        errorbar(4, squeeze(mean_var_area(1,g,v,4)), squeeze(stderr_var_area(1,g,v,4)),'Color',Colors_area{4})
        plot(4, squeeze(mean_var_area(1,g,v,4)),'Linestyle','none','Marker','o','MarkerFaceColor',Colors_area{4},'MarkerEdgeColor',Colors_area{4})
        xticks([1 2 3 4])
        xticklabels({'AUDd','AUDp','AUDv','TEa'})

        box off
        ylabel(v_label(v))
        set(subplot(2,4, v),'YScale',Yscale_all(v));
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
directory = 'Z:\Shared\Benne\Praegel_et_al_2024\praegel_et_al_final\figures';
filename = 'SuppFig412_adolescence'
filepath = fullfile(directory, filename);
recent_figure = figure(1);
saveas(recent_figure, filepath, 'svg');

filename = 'SuppFig412_adult'
filepath = fullfile(directory, filename);
recent_figure = figure(2);
saveas(recent_figure, filepath, 'svg');
%% 
clearvars -except directory GNG_rec_all_cell GNG_analysis 
FR_array = GNG_analysis.frArray_Smoothed_mod ; % only auditory neurons 
bin_array = GNG_analysis.binArray_mod ; % only auditory neurons 
save("FR_array.mat","FR_array")
save("bin_array.mat","bin_array")
