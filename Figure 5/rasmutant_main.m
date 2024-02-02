% Plots included: Fig. 6D, Supplementary Fig. S6, Fig. 6B
%%
%clear all 
clc

% Load baseline parameter values
paramlist;
paramlist_krasmutant; %16-fold decrease in kRhydro
%% Check concentrations for MAPK species
timeSpan = 0:1:60;
%WT Ras w/o sorafenib
[T,~,~,params_minsor,allNames,allValues_wt]             = fullEGFR9_onemodel(timeSpan, yinit, params, 'min_sor', 'no');
%WT Ras w/ sorafenib
[~,~,~,params_plussor,~,~]                              = fullEGFR9_onemodel(timeSpan, yinit, params, 'plus_sor', 'no');
%mutant w/o sorafenib
[~,~,~,params_mutant,~,allValues_mutant]                = fullEGFR9_onemodel(timeSpan, yinit_mutant, params_mutant,'min_sor', 'no');
%Define species of interest
species_list                                            = [1; 2; 3; 4; 5; 6; 7; 8; 9;10; 11; 12; 13; 14;15;16;17;18;19;20;21;22;23;24;...
                                                           25; 26; 27; 28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46]; 

species_names                                           = {'Ras-iBRaf';'Ras-Raf1-iBRaf-tetramer';'Ras-BRaf' ; 'Ras-nfpRaf1' ; 'iRaf1' ; ...
                                                           'ERK'; 'nfpSOS' ; 'Ras-pRaf1' ; 'BRaf-iBRaf-dimer' ; 'Ras-pRaf1-iBRaf_tetramer' ;...
                                                           'Ras-pRaf1-tetramer' ; 'Ras-nfpiRaf1' ; 'Ras-GTP' ; 'nfpiRaf1' ; 'Raf1' ; ...
                                                           'Ras-BRaf-pRaf1-tetramer'; 'Ras-pRaf1-Raf1-tetramer'; 'nfpRaf1' ; 'Ras-nfpBRaf'; ...
                                                           'iBRaf'; 'pMEK'; 'mE'; 'mEL'; 'pRaf1'; 'EG2'; 'Ras-iRaf1'; 'Ras-nfpiBRaf'; 'MEK'; ...
                                                           'Ras-Braf-iRaf1-tetramer'; 'Ras-Raf1'; 'mELmEL'; 'nfpiBRaf'; 'BRaf'; 'SOS'; 'nfpBRaf';...
                                                           'GRB2-SOS'; 'GRB2'; 'Ras-iRaf1-tetramer'; 'iBRaf-dimer'; 'Ras-GDP' ; 'E' ; ...
                                                           'Ras-BRaf-Raf1-tetramer' ; 'BRaf-dimer'; 'EG2SOS' ; 'pERK' ; ...
                                                           'Ras-iRaf1-iBRaf1-tetramer'};
species_group                                           = [species_list string(species_names)];

species_names_cat                                       = categorical(species_names);
species_names_cat                                       = reordercats(species_names_cat,{'Ras-iBRaf';'Ras-Raf1-iBRaf-tetramer';'Ras-BRaf' ; 'Ras-nfpRaf1' ; 'iRaf1' ; ...
                                                           'ERK'; 'nfpSOS' ; 'Ras-pRaf1' ; 'BRaf-iBRaf-dimer' ; 'Ras-pRaf1-iBRaf_tetramer' ;...
                                                           'Ras-pRaf1-tetramer' ; 'Ras-nfpiRaf1' ; 'Ras-GTP' ; 'nfpiRaf1' ; 'Raf1' ; ...
                                                           'Ras-BRaf-pRaf1-tetramer'; 'Ras-pRaf1-Raf1-tetramer'; 'nfpRaf1' ; 'Ras-nfpBRaf'; ...
                                                           'iBRaf'; 'pMEK'; 'mE'; 'mEL'; 'pRaf1'; 'EG2'; 'Ras-iRaf1'; 'Ras-nfpiBRaf'; 'MEK'; ...
                                                           'Ras-Braf-iRaf1-tetramer'; 'Ras-Raf1'; 'mELmEL'; 'nfpiBRaf'; 'BRaf'; 'SOS'; 'nfpBRaf';...
                                                           'GRB2-SOS'; 'GRB2'; 'Ras-iRaf1-tetramer'; 'iBRaf-dimer'; 'Ras-GDP' ; 'E' ; ...
                                                           'Ras-BRaf-Raf1-tetramer' ; 'BRaf-dimer'; 'EG2SOS' ; 'pERK' ; ...
                                                           'Ras-iRaf1-iBRaf1-tetramer'});
%sort based on pathway order
pway_order                                              = {'mE';'mEL';'mELmEL';'E';'EG2';'EG2SOS';...
                                                           'SOS';'GRB2-SOS';'GRB2';'nfpSOS';...
                                                           'Ras-GDP';'Ras-GTP';'BRaf';'BRaf-dimer';'Raf1';'Ras-Raf1';'Ras-BRaf';'Ras-pRaf1';...
                                                           'Ras-nfpBRaf';'Ras-nfpRaf1';'Ras-BRaf-pRaf1-tetramer';'Ras-BRaf-Raf1-tetramer';... 
                                                           'Ras-pRaf1-Raf1-tetramer';'Ras-pRaf1-tetramer';... 
                                                           'pRaf1';'nfpRaf1';'nfpBRaf';...
                                                           'MEK';'pMEK';'ERK';'pERK';...
                                                           'Ras-iRaf1-iBRaf1-tetramer';'nfpiBRaf';'iRaf1' ;...
                                                           'Ras-nfpiRaf1' ;'Ras-iBRaf';'Ras-iRaf1';'iBRaf-dimer'; ...
                                                           'BRaf-iBRaf-dimer' ;'Ras-nfpiBRaf';'Ras-Braf-iRaf1-tetramer';...
                                                           'Ras-Raf1-iBRaf-tetramer';'Ras-iRaf1-tetramer';'nfpiRaf1';...
                                                           'Ras-pRaf1-iBRaf_tetramer';'iBRaf'};

[~,pway_order_num]                                      = ismember(pway_order, species_names);
pway_order_group                                        = [pway_order_num string(pway_order)];

% Check difference in dynamics between Ras mutant and WT Ras
figure
tiledlayout(4,6, 'Padding', 'none', 'TileSpacing', 'compact'); 
for i=1:length(pway_order_num)/2
    %subplot(5,10,i)
    nexttile
    plot(T,allValues_wt(:,pway_order_num(i)), 'LineWidth', 1);
    hold on
    plot(T,allValues_mutant(:,pway_order_num(i)), 'LineWidth', 1);
    species_idx = pway_order(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end
figure
tiledlayout(4,6, 'Padding', 'none', 'TileSpacing', 'compact'); 
for i=length(pway_order_num)/2+1:length(pway_order_num)
    %subplot(5,10,i)
    nexttile
    plot(T,allValues_wt(:,pway_order_num(i)), 'LineWidth', 1);
    hold on
    plot(T,allValues_mutant(:,pway_order_num(i)), 'LineWidth', 1);
    species_idx = pway_order(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end
% checked sorafenib is 0, but iBRaf and iRaf species are not 0
%{
for i=1:size(allValues_mutant,1)
    for j=1:size(allValues_mutant,2)
        if allValues_mutant(i,j) < 0
            allValues_mutant(i,j) = 0;
        end
    end
end
%}
%min-max scaling 
allValues_mutant_norm                                   = zeros(size(allValues_mutant,1), size(allValues_mutant,2)); 
allValues_wt_norm                                       = zeros(size(allValues_wt,1), size(allValues_wt,2)); 
for i=1:length(species_list)
    allValues_mutant_norm(:,i)                          = (allValues_mutant(:,species_list(i)) - min(allValues_mutant(:,species_list(i)))) ./ (max(allValues_mutant(:,species_list(i))) - min(allValues_mutant(:,species_list(i))));
    allValues_wt_norm(:,i)                              = (allValues_wt(:,species_list(i)) - min(allValues_wt(:,species_list(i)))) ./ (max(allValues_wt(:,species_list(i))) - min(allValues_wt(:,species_list(i))));
end

% Check difference in dynamics between Ras mutant and WT Ras
figure
for i=1:length(species_list)
    subplot(5,10,i)
    plot(T,allValues_wt_norm(:,i));
    hold on
    plot(T,allValues_mutant_norm(:,i));
    species_idx = species_names_cat(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end


delta_values = zeros(length(species_list),1);
for i=1:length(species_list)
    %delta_values(i)                                     = trapz(allValues_mutant_norm(:,species_list(i))) ./ trapz(allValues_wt_norm(:,species_list(i)));
    delta_values(i)                                     = trapz(allValues_mutant(:,species_list(i))) ./ trapz(allValues_wt(:,species_list(i)));
end
delta_group                                             = [delta_values string(species_names)];
sorted_delta_group                                      = sortrows(delta_group,1,'descend');
sorted_delta_numbers                                    = str2double(sorted_delta_group(:,1))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
figure
bar(sorted_delta_numbers);
hold on
xlim=get(gca,'xlim');
%plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:length(sorted_delta_numbers),'xTickLabel',sorted_delta_group(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta_numbers),diag(sorted_delta_numbers),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250


pway_delta_values                                       = zeros(length(pway_order_num),1);
for i=1:length(pway_order)
    %pway_delta_values(i)                                = trapz(allValues_mutant_norm(:,pway_order_num(i))) ./ trapz(allValues_wt_norm(:,pway_order_num(i)));
    %pway_delta_values(i)                                = trapz(allValues_mutant(:,pway_order_num(i))) ./ trapz(allValues_wt(:,pway_order_num(i)));
    pway_delta_values(i)                                = abs(trapz(allValues_mutant_norm(:,pway_order_num(i))) - trapz(allValues_wt_norm(:,pway_order_num(i))));

end

%high_delta_index                                        = find(pway_delta_values > 1.5);   
high_delta_index                                        = find(pway_delta_values > 10);   

[delta_group_pway_nonzero_idx, ~]                       = find(pway_delta_values > 0);  

figure
bar(pway_delta_values(delta_group_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
%plot(xlim,[1.5 1.5],'-.k','LineWidth',1);
plot(xlim,[10 10],'-.k','LineWidth',1);
%plot(xlim,[0.5 0.5],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
%ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_pway_nonzero_idx),'xTickLabel',pway_order_group(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_pway_nonzero_idx),diag(pway_delta_values(delta_group_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end
%ylim([0 20]);
set(gcf,'Position',[100 100 750 500]);
%saveas(gcf,[pwd '/Plots/timeint_changes.pdf']);


figure
bar(pway_delta_values(delta_group_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
%plot(xlim,[1.5 1.5],'-.k','LineWidth',1);
plot(xlim,[10 10],'-.k','LineWidth',1);
%plot(xlim,[0.5 0.5],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
%ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_pway_nonzero_idx),'xTickLabel',pway_order_group(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_pway_nonzero_idx),diag(pway_delta_values(delta_group_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [153/255 153/255 255/255]; %227, 111, 71
end
%ylim([0 20]);
set(gcf,'Position',[100 100 750 500]);

%% with 10-fold increase in Raf1 expression:
yinit_raf_highexp                                                   = yinit_mutant;
%yinit_raf_highexp(15)                                               = yinit_mutant(15)*10; % RAF1 increase 10-fold
%yinit_raf_highexp(33)                                               = yinit_mutant(33)*10; % BRAF increase 10-fold
yinit_raf_highexp(34)                                               = yinit_mutant(34)*10; % SOS increase 10-fold

params_raf_highexp                                                  = params;
%params_raf_highexp(115)                                             = params_raf_highexp(115) * 10; % RAF1 increase 10-fold
%params_raf_highexp(131)                                            = params_raf_highexp(131) * 10;
params_raf_highexp(11)                                              = params_raf_highexp(11) * 10; % SOS increase 10-fold


params_mutant_raf_highexp                                           = params_mutant;
%params_mutant_raf_highexp(115)                                      = params_mutant_raf_highexp(115) * 10; % RAF1 increase 10-fold
%params_mutant_raf_highexp(131)                                      = params_mutant_raf_highexp(131) * 10;
params_mutant_raf_highexp(11)                                      = params_mutant_raf_highexp(11) * 10; % SOS increase 10-fold


% WT-Ras dynamics
[~,~,~,out_params_wt_raf_highexp,~,allValues_wt_raf_highexp]        = fullEGFR9_onemodel(timeSpan,yinit_raf_highexp, params_raf_highexp,'min_sor','no'); %wild-type RAS with 10-fold increase in [RAF1]
%[~,~,~,out_params_wt_raf_highexp,~,allValues_wt_raf_highexp]        = fullEGFR9_onemodel(timeSpan,yinit, params,'min_sor','no'); %wild-type RAS with 10-fold increase in [RAF1]

% Mutant Ras dynamics
[~,~,~,out_params_mutant_raf_highexp,~,allValues_mut_raf_highexp]   = fullEGFR9_onemodel(timeSpan,yinit_raf_highexp, params_mutant_raf_highexp,'min_sor','no'); %KRAS G12V with 10-fold increase in [RAF1]
% checked: only hydrolysis rate constant is different between
% wt_raf_highexp and mutant_raf_highexp, init_raf1 12000 * 10

%error tolerance is not working
for i=1:size(allValues_wt_raf_highexp,1)
    for j=1:size(allValues_wt_raf_highexp,2)
        if allValues_wt_raf_highexp(i,j) < 1e-9
            allValues_wt_raf_highexp(i,j) = 0;
        end
    end
end
for i=1:size(allValues_mut_raf_highexp,1)
    for j=1:size(allValues_mut_raf_highexp,2)
        if allValues_mut_raf_highexp(i,j) < 1e-9
            allValues_mut_raf_highexp(i,j) = 0;
        end
    end
end


figure
for i=1:length(species_list)
    subplot(5,10,i)
    plot(T,allValues_wt_raf_highexp(:,species_list(i)));
    hold on
    plot(T,allValues_mut_raf_highexp(:,species_list(i)));
    species_idx = species_names_cat(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end


figure
tiledlayout(4,6, 'Padding', 'none', 'TileSpacing', 'compact'); 
for i=1:length(pway_order_num)/2
    %subplot(5,10,i)
    nexttile
    plot(T,allValues_wt_raf_highexp(:,pway_order_num(i)), 'LineWidth', 1);
    hold on
    plot(T,allValues_mut_raf_highexp(:,pway_order_num(i)), 'LineWidth', 1);
    species_idx = pway_order(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end
figure
tiledlayout(4,6, 'Padding', 'none', 'TileSpacing', 'compact'); 
for i=length(pway_order_num)/2+1:length(pway_order_num)
    %subplot(5,10,i)
    nexttile
    plot(T,allValues_wt_raf_highexp(:,pway_order_num(i)), 'LineWidth', 1);
    hold on
    plot(T,allValues_mut_raf_highexp(:,pway_order_num(i)), 'LineWidth', 1);
    species_idx = pway_order(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end

%min-max scaling 
allValues_rafinc_mutant_norm                                        = zeros(size(allValues_mut_raf_highexp,1), size(allValues_mut_raf_highexp,2)); 
allValues_rafinc_wt_norm                                            = zeros(size(allValues_wt_raf_highexp,1), size(allValues_wt_raf_highexp,2)); 
for i=1:length(species_list)
    allValues_rafinc_mutant_norm(:,i)                               = (allValues_mut_raf_highexp(:,species_list(i)) - min(allValues_mut_raf_highexp(:,species_list(i)))) ./ (max(allValues_mut_raf_highexp(:,species_list(i))) - min(allValues_mut_raf_highexp(:,species_list(i))));
    allValues_rafinc_wt_norm(:,i)                                   = (allValues_wt_raf_highexp(:,species_list(i)) - min(allValues_wt_raf_highexp(:,species_list(i)))) ./ (max(allValues_wt_raf_highexp(:,species_list(i))) - min(allValues_wt_raf_highexp(:,species_list(i))));
end

% Check difference in dynamics between Ras mutant and WT Ras
figure
for i=1:length(species_list)
    subplot(5,10,i)
    plot(T,allValues_rafinc_wt_norm(:,species_list(i)));
    hold on
    plot(T,allValues_rafinc_mutant_norm(:,species_list(i)));
    species_idx = species_names_cat(i);
    title(species_idx);
    if i==1
        AX=legend('WT Ras','Mutant Ras','location','northeast');
    end
end
% Plot in descending order
delta_values_rafinc                                                 = zeros(length(species_list),1);
for i=1:length(species_list)
    %delta_values_rafinc(i)                                          = trapz(allValues_rafinc_mutant_norm(:,species_list(i))) ./ trapz(allValues_rafinc_wt_norm(:,species_list(i)));
    %delta_values_rafinc(i)                                          = trapz(allValues_mut_raf_highexp(:,species_list(i))) ./ trapz(allValues_wt_raf_highexp(:,species_list(i)));
    delta_values_rafinc(i)                                          = abs(trapz(allValues_rafinc_mutant_norm(:,species_list(i))) - trapz(allValues_rafinc_wt_norm(:,species_list(i))));

end
delta_group_rafinc                                                  = [delta_values_rafinc string(species_names)];
sorted_delta_group_rafinc                                           = sortrows(delta_group_rafinc,1,'descend');
sorted_delta_numbers_rafinc                                         = str2double(sorted_delta_group_rafinc(:,1));
% Highlight everything with above 1 fold change
high_delta_index                                                    = find(delta_values_rafinc > 1);

figure
bar(sorted_delta_numbers_rafinc);
hold on
xlim=get(gca,'xlim');
xlabel('Node','FontWeight','bold');
ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:length(sorted_delta_numbers_rafinc),'xTickLabel',sorted_delta_group_rafinc(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta_numbers_rafinc),diag(sorted_delta_numbers_rafinc),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end

% Plot in pathway order
delta_values_rafinc_pway                                                 = zeros(length(pway_order_num),1);
for i=1:length(pway_order_num)
    %delta_values_rafinc_pway(i)                                          = trapz(allValues_rafinc_mutant_norm(:,pway_order_num(i))) ./ trapz(allValues_rafinc_wt_norm(:,pway_order_num(i)));
    %delta_values_rafinc_pway(i)                                          = trapz(allValues_mut_raf_highexp(:,pway_order_num(i))) ./ trapz(allValues_wt_raf_highexp(:,pway_order_num(i)));
     delta_values_rafinc_pway(i)                                          = abs(trapz(allValues_rafinc_mutant_norm(:,pway_order_num(i))) - trapz(allValues_rafinc_wt_norm(:,pway_order_num(i))));

end
delta_group_rafinc_pway                                                  = [delta_values_rafinc_pway string(pway_order_num)];
% Highlight everything with above 1.5 fold change
%high_delta_index                                                    = find(delta_values_rafinc_pway > 1.5);    
high_delta_index                                                    = find(delta_values_rafinc_pway > 10);    

[delta_group_rafinc_pway_nonzero_idx, vals]                         = find(delta_values_rafinc_pway > 0);  

figure
bar(delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
plot(xlim,[10 10],'-.k','LineWidth',1);
%plot(xlim,[1.5 1.5],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
%ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_rafinc_pway_nonzero_idx),'xTickLabel',pway_order_group(delta_group_rafinc_pway_nonzero_idx,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_rafinc_pway_nonzero_idx),diag(delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end
%ylim([0 20]);
set(gcf,'Position',[100 100 750 500]);
saveas(gcf,[pwd '/Plots/rafinc_timeint_changes.pdf']);

delta_baseline_inc = (delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx)- pway_delta_values(delta_group_pway_nonzero_idx));
high_delta_index                                                    = find(delta_baseline_inc > 5);    
%[delta_group_rafinc_pway_nonzero_idx, vals]                         = find(delta_baseline_inc > 5);  

figure
bar(delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
plot(xlim,[10 10],'-.k','LineWidth',1);
%plot(xlim,[1.5 1.5],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
%ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_rafinc_pway_nonzero_idx),'xTickLabel',pway_order_group(delta_group_rafinc_pway_nonzero_idx,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_rafinc_pway_nonzero_idx),diag(delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end
%ylim([0 20]);
set(gcf,'Position',[100 100 750 500]);

%% Previous versions
[delta_Sorted, delta_index]                                         = sort(delta_species,'descend');
sorted_delta                                                        = delta_species(delta_index);
high_delta                                                          = sorted_delta(sorted_delta > 2);
high_delta_index                                                    = find(sorted_delta > 2);
%remove NaN values from the difference 
[nan_row, nan_col]                                                  = find(isnan(sorted_delta));
sorted_delta_wo_nan                                                 = sorted_delta(~isnan(sorted_delta));
species_names_wo_nan                                                = species_names;
species_names_wo_nan(nan_row) = [];


figure
bar(sorted_delta_wo_nan);
hold on
xlim=get(gca,'xlim');
%plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Time-integrated, normalized concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:length(sorted_delta_wo_nan),'xTickLabel',species_names_wo_nan,'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta_wo_nan),diag(sorted_delta_wo_nan),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    %bh(k).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end




nodelabel = {'GTP-RAS'; 'Total GTP-bound RAS'; ...
    'Membrane-bound RAF'; ...
    'BRaf-dimer'; ...
    'Ras-BRaf-Raf1-tetramer'; ...
    'Ras-BRaf-pRaf1-tetramer'; ...
    'Ras-pRaf1'; ...
    'Ras-pRaf1-Raf1-tetramer'; ...
    'Ras-pRaf1-tetramer';...
    'pRaf1'; ...
    'pMEK'; ...
    'pERK'; ...
    'Ras-Raf1'};
nodelabel = categorical(nodelabel);
nodelabel2 = reordercats(nodelabel,{'GTP-RAS'; 'Total GTP-bound RAS'; ...
    'Membrane-bound RAF'; ...
    'BRaf-dimer'; ...
    'Ras-BRaf-Raf1-tetramer'; ...
    'Ras-BRaf-pRaf1-tetramer'; ...
    'Ras-pRaf1'; ...
    'Ras-pRaf1-Raf1-tetramer'; ...
    'Ras-pRaf1-tetramer';...
    'pRaf1'; ...
    'pMEK'; ...
    'pERK'; ...
    'Ras-Raf1'});
%{
figure
bar(nodelabel2,delta_species);
hold on 
xlim=get(gca,'xlim');plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node');
ylabel({'Fold change in maximum outputs';'Mutant / WT'});
title('Effect of 16-fold decrease in RAS-GTP hydrolysis rate');
%}
%% Fig. 6C
[delta_Sorted, delta_index] = sort(delta_species,'descend');
sorted_delta = delta_species(delta_index);
high_delta = sorted_delta(sorted_delta > 2);
high_delta_index = find(sorted_delta > 2);
figure
bar(sorted_delta);
hold on
xlim=get(gca,'xlim');
plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Maximum concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:13,'xTickLabel',nodelabel2(delta_index),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta),diag(sorted_delta),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(k).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end





%% Test the effect of increasing RAF1 expression in wild-type RAS & KRAS mutant (Fig. 6D)
yinit_raf_highexp = yinit_mutant;
yinit_raf_highexp(15) = yinit_mutant(15)*10; % RAF1 increase 10-fold

[~,~,~,params_raf_highexp,~,allValues_raf_highexp] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_mutant,'min_sor', 'no'); %KRAS G12V with 10-fold increase in [RAF1]
[~,~,~,params_wt_raf_highexp,~,allValues_wt_raf_highexp] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params,'min_sor', 'no'); %wild-type RAS with 10-fold increase in [RAF1]

h1 = figure;
set(h1,'Position',[100 100 470 300]);
%plot(T,allValues_wt(:,45),'LineWidth', 1.5,'Color',[0 0 0]); %WT
plot(T,allValues_wt(:,45),'LineWidth', 1.5,'Color',[0/255 191/255 196/255]); %WT
hold on
plot(T,allValues_wt_raf_highexp(:,45),'--','LineWidth', 1.5,'Color',[0/255 191/255 196/255]); %WT w/ RAF increase
%plot(T,allValues_mutant(:,45),'LineWidth', 1.5, 'Color', [255/255 0 0]); %KRAS mutant 
plot(T,allValues_mutant(:,45),'LineWidth', 1.5, 'Color', [248/255 118/255 109/255]); %KRAS mutant 
%plot(T,allValues_raf_highexp(:,45), '--', 'LineWidth', 1.5,'Color', [255/255 0 0]); %KRAS mutant w/ RAF increase
plot(T,allValues_raf_highexp(:,45), '--', 'LineWidth', 1.5,'Color', [248/255 118/255 109/255]); %KRAS mutant w/ RAF increase
hold off
legend('Wild-type RAS', 'Wild-type RAS, 10-fold increase in [RAF1]', 'KRAS G12V mutant', 'KRAS G12V mutant, 10-fold increase in [RAF1]','Location','northoutside','FontSize',8);
xlabel('Time (min)','FontSize',8, 'FontWeight','bold');
ylabel('pERK (molec/cell)','FontSize',8, 'FontWeight','bold');
ylim([0 3.8*10^5]);
ax = gca;
ax.FontSize = 6; 
set(gcf,'Position',[100 100 750 500])
savefig('mutant_wt_raf_init_change.fig');

fold_diff1 = abs(trapz(T,allValues_wt(:,45)) - trapz(T,allValues_mutant(:,45)))
fold_diff2 = abs(trapz(T,allValues_wt_raf_highexp(:,45)) - trapz(T,allValues_raf_highexp(:,45)))
rafinc_effect = fold_diff2/fold_diff1

%% Test feedback role (nfpSOS, nfpBR1)  
params_mutant_nofdback     = params_mutant;
params_mutant_nofdback(10) = 0; % param(10) is 'knfpSOS'
params_nofdback            = params;
params_nofdback(10)        = 0;
[~,~,~,params_raf_highexp_nofdbk,~,allValues_raf_highexp_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_mutant_nofdback,'min_sor', 'no'); %KRAS G12V with 10-fold increase in [RAF1]
[~,~,~,params_wt_raf_highexp_nofdbk,~,allValues_wt_raf_highexp_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_nofdback,'min_sor', 'no'); %wild-type RAS with 10-fold increase in [RAF1]

h1 = figure;
set(h1,'Position',[100 100 380 300]);
%plot(T,allValues_wt(:,45),'LineWidth', 1.5,'Color',[0 0 0]); %WT
plot(T,allValues_wt(:,45),'LineWidth', 1.5,'Color',[0/255 191/255 196/255]); %WT
hold on
plot(T,allValues_wt_raf_highexp_nofdbk(:,45),'--','LineWidth', 1.5,'Color',[0/255 191/255 196/255]); %WT w/ RAF increase
plot(T,allValues_mutant(:,45),'LineWidth', 1.5, 'Color', [248/255 118/255 109/255]); %KRAS mutant 
plot(T,allValues_raf_highexp_nofdbk(:,45), '--', 'LineWidth', 1.5,'Color', [248/255 118/255 109/255]); %KRAS mutant w/ RAF increase
hold off
legend('Wild-type RAS', 'Wild-type RAS, 10-fold increase in [RAF1], no fdback', 'KRAS G12V mutant', 'KRAS G12V mutant, 10-fold increase in [RAF1], no fdback','Location','northoutside','FontSize',8);
xlabel('Time (min)','FontSize',12, 'FontWeight','bold');
ylabel('pERK (molec/cell)','FontSize',12, 'FontWeight','bold');
ylim([0 4*10^5]);
ax = gca;
ax.FontSize = 8; 
set(gcf,'Position',[100 100 750 500])
savefig('mutant_wt_raf_init_change.fig');

fold_diff1_nofdbck    = abs(trapz(T,allValues_wt(:,45)) - trapz(T,allValues_mutant(:,45)))
fold_diff2_nofdbck    = abs(trapz(T,allValues_wt_raf_highexp_nofdbk(:,45)) - trapz(T,allValues_raf_highexp_nofdbk(:,45)))
rafinc_effect_nofdbck = fold_diff2_nofdbck/fold_diff1_nofdbck

%% time integrated pERK diff. between WT vs. Mutant
bar_x = categorical({'12,000', '120,000', '120,000 SOS off'});
bar_x = reordercats(bar_x, {'12,000', '120,000', '120,000 SOS off'});
bar_y = [fold_diff1 fold_diff2 fold_diff2_nofdbck];
figure
set(gcf,'Position',[100 100 300 400])
bar(bar_x, bar_y, 'BarWidth',0.6);
ylim([0 max(bar_y)*1.2]);

%% Test the effect of Kd (Supplementary Fig. S6)
new_kd = [1.251011e+04*0.01 1.251011e+04*0.1 1.251011e+04];
new_params_EGFRdimer_minsor = zeros(length(T), 1);
new_params_pEGFR_minsor = zeros(length(T), 1);
new_params_RasGTP_minsor = zeros(length(T), 1);
new_params_membRAF1_minsor = zeros(length(T), 1);
new_params_pMEK_minsor = zeros(length(T), 1);
new_params_pERK_minsor = zeros(length(T), 1);

for i=1:length(new_kd)
    new_params = params;
    new_params(2) = new_kd(i);
    [T,~,~,new_params_minsor,~,new_params_allValues_minsor]   = fullEGFR9_onemodel(timeSpan, yinit, new_params, 'min_sor', 'no');
    [~,~,~,new_params_plussor,~,new_params_allValues_plussor] = fullEGFR9_onemodel(timeSpan, yinit, new_params, 'plus_sor', 'no');
    
    
    new_params_EGFRdimer_minsor(:,i)                          = new_params_allValues_minsor(:,31);
    new_params_pEGFR_minsor(:,i)                              = new_params_allValues_minsor(:,41);
    new_params_RasGTP_minsor(:,i)                             = new_params_allValues_minsor(:,184);
    new_params_membRAF1_minsor(:,i)                           = new_params_allValues_minsor(:,157);
    new_params_pMEK_minsor(:,i) = new_params_allValues_minsor(:,21);
    new_params_pERK_minsor(:,i) = new_params_allValues_minsor(:,45);
end
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
for i=1:length(new_kd)
    set(0, 'CurrentFigure', f1)
    plot(T, new_params_EGFRdimer_minsor(:,i),'LineWidth',1.2)
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF');
    title('EGFR dimer');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('EGFR dimer (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f2)
    plot(T, new_params_pEGFR_minsor(:,i),'LineWidth',1.2)
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF','Location','Southeast');
    title('Phosphorylated EGFR');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pEGFR (molec/cell)','FontWeight','bold');
    hold on

    set(0, 'CurrentFigure', f3)
    plot(T, new_params_RasGTP_minsor(:,i),'LineWidth',1.2)
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF');
    title('GTP-bound RAS');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('GTP-bound RAS (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f4)
    plot(T, new_params_membRAF1_minsor(:,i),'LineWidth',1.2)
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF');
    title('Membrane-bound RAF1');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('Membrane-bound RAF1 (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f5)
    plot(T, new_params_pMEK_minsor(:,i),'LineWidth',1.2)
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF');
    title('pMEK');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pMEK (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f6)
    plot(T, new_params_pERK_minsor(:,i),'LineWidth',1.2)
    title('pERK');
    legend('100-fold lower affinity ligand than EGF', '10-fold lower affinity ligand than EGF', 'EGF');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pERK (molec/cell)','FontWeight','bold');
    hold on
end

%% Test the effect of changing RAF1 expression (Supplementary Fig. S6)
new_init_param = [12000.0*100 12000.0*10 12000.0];
new_init_params_EGFRdimer_minsor = zeros(length(T), 1);
new_init_params_pEGFR_minsor = zeros(length(T), 1);
new_init_params_RasGTP_minsor = zeros(length(T), 1);
new_init_params_membRAF1_minsor = zeros(length(T), 1);
new_init_params_pRAF1_minsor = zeros(length(T), 1);
new_init_params_pMEK_minsor = zeros(length(T), 1);
new_init_params_pERK_minsor = zeros(length(T), 1);
RAF1diff = cell(2,1);

for i=1:length(new_init_param)
    new_init = yinit;
    new_init(15) = new_init_param(i);
    new_init_params = params;
    new_init_params(141) = new_init_param(i);
    
    [T,~,~,new_params_minsor,~,new_init_params_allValues_minsor] = fullEGFR9_onemodel(timeSpan, new_init, new_init_params, 'min_sor', 'no');
    [~,~,~,new_params_plussor,~,new_init_params_allValues_plussor] = fullEGFR9_onemodel(timeSpan, new_init, new_init_params, 'plus_sor', 'no');
    
    RAF1diff{i} = new_init_params_allValues_minsor-allValues_wt;
    
    new_init_params_EGFRdimer_minsor(:,i) = new_init_params_allValues_minsor(:,31);
    new_init_params_pEGFR_minsor(:,i) = new_init_params_allValues_minsor(:,41);
    new_init_params_RasGTP_minsor(:,i) = new_init_params_allValues_minsor(:,184);
    new_init_params_membRAF1_minsor(:,i) = new_init_params_allValues_minsor(:,157);
    new_init_params_pRAF1_minsor(:,i) = new_init_params_allValues_minsor(:,24);
    new_init_params_pMEK_minsor(:,i) = new_init_params_allValues_minsor(:,21);
    new_init_params_pERK_minsor(:,i) = new_init_params_allValues_minsor(:,45);
end


f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
f7 = figure;
C = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]}; 
for i=1:length(new_init_param)
    set(0, 'CurrentFigure', f1)
    plot(T, new_init_params_EGFRdimer_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell');
    title('EGFR dimer');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('EGFR dimer (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f2)
    plot(T, new_init_params_pEGFR_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell');
    title('Phosphorylated EGFR');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pEGFR (molec/cell)','FontWeight','bold');
    hold on

    set(0, 'CurrentFigure', f3)
    plot(T, new_init_params_RasGTP_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell');
    title('GTP-bound RAS');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('GTP-bound RAS (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f4)
    plot(T, new_init_params_membRAF1_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell');
    title('Membrane-bound RAF1');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('Membrane-bound RAF1 (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f5)
    plot(T, new_init_params_pRAF1_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell');
    title('pRAF1');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pRAF1 (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f6)
    plot(T, new_init_params_pMEK_minsor(:,i),'LineWidth',1.2, 'Color',C{i})
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell','Location','southeast');
    title('pMEK');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pMEK (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f7)
    plot(T, new_init_params_pERK_minsor(:,i),'LineWidth',1.2,'Color',C{i})
    title('pERK');
    legend('100-fold decrease in [RAF1]','10-fold decrease in [RAF1]','[RAF1] = 12,000 molec/cell','Location','southeast');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pERK (molec/cell)','FontWeight','bold');
    ylim([0 max(new_init_params_pERK_minsor,[],'all')+max(new_init_params_pERK_minsor,[],'all')*0.1])
    hold on
end
%saveas(f6,[pwd '/Plots/pmek_raf_init_change.fig']);
%saveas(f6,[pwd '/Plots/pmek_raf_init_change.pdf']);
%saveas(f7,[pwd '/Plots/perk_raf_init_change.fig']);
%saveas(f7,[pwd '/Plots/perk_raf_init_change.pdf']);

%% Test the effect of low EGFR expression (not included)
new_init_param = [93000.0*0.01 93000.0];
new_init_params_EGFRdimer_minsor = zeros(length(T), 1);
new_init_params_pEGFR_minsor = zeros(length(T), 1);
new_init_params_RasGTP_minsor = zeros(length(T), 1);
new_init_params_membRAF1_minsor = zeros(length(T), 1);
new_init_params_pRAF1_minsor = zeros(length(T), 1);
new_init_params_pMEK_minsor = zeros(length(T), 1);
new_init_params_pERK_minsor = zeros(length(T), 1);
EGFRdiff = cell(2,1);

for i=1:length(new_init_param)
    new_init = yinit;
    new_init(22) = new_init_param(i);
    new_init_params = params;
    new_init_params(141) = new_init_param(i);
    
    [T,~,~,new_params_minsor,~,new_init_params_allValues_minsor] = fullEGFR9_onemodel(timeSpan, new_init, new_init_params, 'min_sor', 'no');
    [~,~,~,new_params_plussor,~,new_init_params_allValues_plussor] = fullEGFR9_onemodel(timeSpan, new_init, new_init_params, 'plus_sor', 'no');
    
    EGFRdiff{i} = new_init_params_allValues_minsor-allValues_wt;
    
    new_init_params_EGFRdimer_minsor(:,i) = new_init_params_allValues_minsor(:,31);
    new_init_params_pEGFR_minsor(:,i) = new_init_params_allValues_minsor(:,41);
    new_init_params_RasGTP_minsor(:,i) = new_init_params_allValues_minsor(:,184);
    new_init_params_membRAF1_minsor(:,i) = new_init_params_allValues_minsor(:,157);
    new_init_params_pRAF1_minsor(:,i) = new_init_params_allValues_minsor(:,24);
    new_init_params_pMEK_minsor(:,i) = new_init_params_allValues_minsor(:,21);
    new_init_params_pERK_minsor(:,i) = new_init_params_allValues_minsor(:,45);
end

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
f7 = figure;
for i=1:length(new_init_param)
    set(0, 'CurrentFigure', f1)
    plot(T, new_init_params_EGFRdimer_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('EGFR dimer');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('EGFR dimer (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f2)
    plot(T, new_init_params_pEGFR_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('Phosphorylated EGFR');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pEGFR (molec/cell)','FontWeight','bold');
    hold on

    set(0, 'CurrentFigure', f3)
    plot(T, new_init_params_RasGTP_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('GTP-bound RAS');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('GTP-bound RAS (molec/cell)','FontWeight','bold');
    hold on

    
    set(0, 'CurrentFigure', f4)
    plot(T, new_init_params_membRAF1_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('Membrane-bound RAF1');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('Membrane-bound RAF1 (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f5)
    plot(T, new_init_params_pRAF1_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('pRAF1');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pRAF1 (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f6)
    plot(T, new_init_params_pMEK_minsor(:,i),'LineWidth',1.2)
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    title('pMEK');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pMEK (molec/cell)','FontWeight','bold');
    hold on
    
    set(0, 'CurrentFigure', f7)
    plot(T, new_init_params_pERK_minsor(:,i),'LineWidth',1.2)
    title('pERK');
    legend('100-fold decrease in [EGFR]','[EGF] = 93000 molec/cell');
    xlabel('Time (min)','FontWeight','bold');
    ylabel('pERK (molec/cell)','FontWeight','bold');
    hold on
end

%Find species that cause pERK to be as sustained as in [EGFR] = 93000 molec/cell
EGFRavgdiff = abs(mean(EGFRdiff{1})); 
[diffval, diffindex] = sort(EGFRavgdiff(:),'descend');
leastdiff = diffindex(2:20);
speciesnames = categorical(allNames(leastdiff));
f1 = figure
bar(speciesnames,diffval(2:20));
%set(gca, 'XtickLabel', speciesnames);
%f1.Position = [10 10 1000 1000];
xtickangle(45);
%%
%[T,~,~,params_minsor,~,allValues_wt] = fullEGFR9_onemodel(timeSpan, yinit, params, 'min_sor');
%[~,~,~,params_plussor,~,~] = fullEGFR9_onemodel(timeSpan, yinit, params, 'plus_sor');
%[T_mutant,~,~,params_mutant,~,allValues_mutant] = fullEGFR9_onemodel(timeSpan, yinit_mutant, params_mutant,'min_sor');

%normalize allValues to their max
allValues_wt_norm = zeros(size(allValues_wt,1),size(allValues_wt,2));
allValues_mutant_norm = zeros(size(allValues_mutant,1),size(allValues_mutant,2));
for i=1:size(allValues_wt,2)
    allValues_wt_norm(:,i) = allValues_wt(:,i)/max(allValues_wt(:,i));
    allValues_mutant_norm(:,i) = allValues_mutant(:,i)/max(allValues_mutant(:,i));
end


%% Identify nodes/species at which RAS hydrolysis increase has the greatest effect

nodes = [%40; %Ras-GDP
    13; %Ras-GTP
    184; %GTP bound RAS
    157; %Membrane bound RAF
    43; %BRaf_dimer
    42; %Ras_BRaf_Raf1_tetramer 
    16;%Ras_BRaf_pRaf1_tetramer 
    8; %Ras_pRaf1 
    17;%Ras_pRaf1_Raf1_tetramer 
    11; %Ras_pRaf1_tetramer
    24; %pRaf1
    21; %pMEK
    45; %pERK 
    30]; %Ras_Raf1
    %132]; %active_Raf]; 

delta_species = zeros(length(nodes),1);
for i=1:length(nodes)
    delta_species(i) = max(allValues_mutant(:,nodes(i))) ./ max(allValues_wt(:,nodes(i)));
end

nodelabel = {'GTP-RAS'; 'Total GTP-bound RAS'; ...
    'Membrane-bound RAF'; ...
    'BRaf-dimer'; ...
    'Ras-BRaf-Raf1-tetramer'; ...
    'Ras-BRaf-pRaf1-tetramer'; ...
    'Ras-pRaf1'; ...
    'Ras-pRaf1-Raf1-tetramer'; ...
    'Ras-pRaf1-tetramer';...
    'pRaf1'; ...
    'pMEK'; ...
    'pERK'; ...
    'Ras-Raf1'};
nodelabel = categorical(nodelabel);
nodelabel2 = reordercats(nodelabel,{'GTP-RAS'; 'Total GTP-bound RAS'; ...
    'Membrane-bound RAF'; ...
    'BRaf-dimer'; ...
    'Ras-BRaf-Raf1-tetramer'; ...
    'Ras-BRaf-pRaf1-tetramer'; ...
    'Ras-pRaf1'; ...
    'Ras-pRaf1-Raf1-tetramer'; ...
    'Ras-pRaf1-tetramer';...
    'pRaf1'; ...
    'pMEK'; ...
    'pERK'; ...
    'Ras-Raf1'});
%{
figure
bar(nodelabel2,delta_species);
hold on 
xlim=get(gca,'xlim');plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node');
ylabel({'Fold change in maximum outputs';'Mutant / WT'});
title('Effect of 16-fold decrease in RAS-GTP hydrolysis rate');
%}
%% Fig. 6C
[delta_Sorted, delta_index] = sort(delta_species,'descend');
sorted_delta = delta_species(delta_index);
high_delta = sorted_delta(sorted_delta > 2);
high_delta_index = find(sorted_delta > 2);
figure
bar(sorted_delta);
hold on
xlim=get(gca,'xlim');
plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Maximum concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:13,'xTickLabel',nodelabel2(delta_index),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta),diag(sorted_delta),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(k).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end
savefig('mutant_wt_foldchange.fig');


nodes2 = [1;
    2;
    3; %GTP bound RAS
    4; %Membrane bound RAF
    5; %BRaf_dimer
    6; %Ras_BRaf_Raf1_tetramer 
    7;%Ras_BRaf_pRaf1_tetramer 
    8; %Ras_pRaf1 
    9;%Ras_pRaf1_Raf1_tetramer 
    10; %Ras_pRaf1_tetramer
    11; %pRaf1
    12; %pMEK
    13; %pERK 
    14;
    15;
    16;
    17;
    18;
    19;
    20;
    21;
    22;
    23;
    24;
    25;
    26;
    27;
    28;
    29;
    30;
    31;
    32;
    33;
    34;
    35;
    36;
    37;
    38;
    39;
    40;
    41;
    42;
    43;
    44;
    45;
    46]; 

nodelabel2 = {'Ras-iBRaf';
	'Ras-Raf1-iBRaf-tetramer';
	'Ras-BRaf' ;
	'Ras-nfpRaf1' ;
	'iRaf1' ;
	'ERK';
	'nfpSOS' ;
	'Ras-pRaf1' ;
	'BRaf-iBRaf-dimer' ;
	'Ras-pRaf1-iBRaf_tetramer' ;
	'Ras-pRaf1-tetramer' ;
	'Ras-nfpiRaf1' ;
	'Ras-GTP' ;
	'nfpiRaf1' ;
	'Raf1' ;
	'Ras-BRaf-pRaf1-tetramer'; 
	'Ras-pRaf1-Raf1-tetramer';
	'nfpRaf1' ;
	'Ras-nfpBRaf';
	'iBRaf';
	'pMEK';
	'mE';
	'mEL';
	'pRaf1';
	'EG2';
	'Ras-iRaf1';
	'Ras-nfpiBRaf';
	'MEK';
	'Ras-Braf-iRaf1-tetramer';
	'Ras-Raf1';
	'mELmEL';
	'nfpiBRaf';
	'BRaf';
	'SOS';
	'nfpBRaf';
	'GRB2-SOS';
	'GRB2';
	'Ras-iRaf1-tetramer';
	'iBRaf-dimer';
	'Ras-GDP' ;
	'E' ;
	'Ras-BRaf-Raf1-tetramer' ;
	'BRaf-dimer';
	'EG2SOS' ;
	'pERK' ;
	'Ras-iRaf1-iBRaf1-tetramer'};

 nodelabel_cat = categorical(nodelabel2);
 nodelabel_cat2 = reordercats(nodelabel_cat,{'Ras-iBRaf';
	'Ras-Raf1-iBRaf-tetramer';
	'Ras-BRaf' ;
	'Ras-nfpRaf1' ;
	'iRaf1' ;
	'ERK';
	'nfpSOS' ;
	'Ras-pRaf1' ;
	'BRaf-iBRaf-dimer' ;
	'Ras-pRaf1-iBRaf_tetramer' ;
	'Ras-pRaf1-tetramer' ;
	'Ras-nfpiRaf1' ;
	'Ras-GTP' ;
	'nfpiRaf1' ;
	'Raf1' ;
	'Ras-BRaf-pRaf1-tetramer'; 
	'Ras-pRaf1-Raf1-tetramer';
	'nfpRaf1' ;
	'Ras-nfpBRaf';
	'iBRaf';
	'pMEK';
	'mE';
	'mEL';
	'pRaf1';
	'EG2';
	'Ras-iRaf1';
	'Ras-nfpiBRaf';
	'MEK';
	'Ras-Braf-iRaf1-tetramer';
	'Ras-Raf1';
	'mELmEL';
	'nfpiBRaf';
	'BRaf';
	'SOS';
	'nfpBRaf';
	'GRB2-SOS';
	'GRB2';
	'Ras-iRaf1-tetramer';
	'iBRaf-dimer';
	'Ras-GDP' ;
	'E' ;
	'Ras-BRaf-Raf1-tetramer' ;
	'BRaf-dimer';
	'EG2SOS' ;
	'pERK' ;
	'Ras-iRaf1-iBRaf1-tetramer'});

nodelabel_minsor = {'Ras-BRaf' ;
	'Ras-nfpRaf1' ;
	'ERK';
	'nfpSOS' ;
	'Ras-pRaf1' ;
	'Ras-pRaf1-tetramer' ;
	'Ras-GTP' ;
	'Raf1' ;
	'Ras-BRaf-pRaf1-tetramer'; 
	'Ras-pRaf1-Raf1-tetramer';
	'nfpRaf1' ;
	'Ras-nfpBRaf';
	'pMEK';
	'mE';
	'mEL';
	'pRaf1';
	'EG2';
	'MEK';
	'Ras-Raf1';
	'mELmEL';
	'BRaf';
	'SOS';
	'nfpBRaf';
	'GRB2-SOS';
	'GRB2';
	'Ras-GDP' ;
	'E' ;
	'Ras-BRaf-Raf1-tetramer' ;
	'BRaf-dimer';
	'EG2SOS' ;
	'pERK'};

 nodelabel_cat_minsor = categorical(nodelabel_minsor);
 nodelabel_cat_minsor = reordercats(nodelabel_cat_minsor,{'Ras-BRaf' ;
	'Ras-nfpRaf1' ;
	'ERK';
	'nfpSOS' ;
	'Ras-pRaf1' ;
	'Ras-pRaf1-tetramer' ;
	'Ras-GTP' ;
	'Raf1' ;
	'Ras-BRaf-pRaf1-tetramer'; 
	'Ras-pRaf1-Raf1-tetramer';
	'nfpRaf1' ;
	'Ras-nfpBRaf';
	'pMEK';
	'mE';
	'mEL';
	'pRaf1';
	'EG2';
	'MEK';
	'Ras-Raf1';
	'mELmEL';
	'BRaf';
	'SOS';
	'nfpBRaf';
	'GRB2-SOS';
	'GRB2';
	'Ras-GDP' ;
	'E' ;
	'Ras-BRaf-Raf1-tetramer' ;
	'BRaf-dimer';
	'EG2SOS' ;
	'pERK'});
nodes2_minsor = [3; %GTP bound RAS
    4; %Membrane bound RAF
    6; %Ras_BRaf_Raf1_tetramer 
    7;%Ras_BRaf_pRaf1_tetramer 
    8; %Ras_pRaf1 
    11; %pRaf1
    13; %pERK 
    15;
    16;
    17;
    18;
    19;
    21;
    22;
    23;
    24;
    25;
    28;
    30;
    31;
    33;
    34;
    35;
    36;
    37;
    40;
    41;
    42;
    43;
    44;
    45]; 


 nodelabel_cat = categorical(nodelabel2);
 nodelabel_cat2 = reordercats(nodelabel_cat,{'Ras-BRaf' ;
	'Ras-nfpRaf1' ;
	'ERK';
	'nfpSOS' ;
	'Ras-pRaf1' ;
	'Ras-pRaf1-tetramer' ;
	'Ras-GTP' ;
	'Raf1' ;
	'Ras-BRaf-pRaf1-tetramer'; 
	'Ras-pRaf1-Raf1-tetramer';
	'nfpRaf1' ;
	'Ras-nfpBRaf';
	'pMEK';
	'mE';
	'mEL';
	'pRaf1';
	'EG2';
	'MEK';
	'Ras-Raf1';
	'mELmEL';
	'BRaf';
	'SOS';
	'nfpBRaf';
	'GRB2-SOS';
	'GRB2';
	'Ras-GDP' ;
	'E' ;
	'Ras-BRaf-Raf1-tetramer' ;
	'BRaf-dimer';
	'EG2SOS' ;
	'pERK'});

delta_species2 = zeros(length(nodes2_minsor),1);
for i=1:length(nodes2_minsor)
    delta_species2(i) = trapz(allValues_mutant(:,nodes2_minsor(i))) ./ max(allValues_wt(:,nodes2_minsor(i)));
end

[delta_Sorted2, delta_index2] = sort(delta_species2,'descend');
sorted_delta2 = delta_species2(delta_index2);
sorted_delta3 = sorted_delta2(~isnan(sorted_delta2))

high_delta2 = sorted_delta3(sorted_delta3 > 2);
high_delta_index2 = find(sorted_delta3 > 2);
figure
bar(sorted_delta3);
hold on
xlim=get(gca,'xlim');
plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Maximum concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
set(gca, 'xTick', 1:46,'xTickLabel',nodelabel_minsor(delta_index2),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta3),diag(sorted_delta3),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index2)
    bh(k).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end

%% Fig. 6C (increase in RAF1 expression effect) for all species
%normalize allValues to their max
allValues_wt_norm           = zeros(size(allValues_wt_raf_highexp,1),size(allValues_wt_raf_highexp,2));
allValues_mutant_norm = zeros(size(allValues_raf_highexp,1),size(allValues_raf_highexp,2));
for i=1:size(allValues_wt_raf_highexp,2)
    allValues_wt_norm(:,i) = allValues_wt_raf_highexp(:,i)./max(allValues_wt_raf_highexp(:,i));
    allValues_mutant_norm(:,i) = allValues_raf_highexp(:,i)./max(allValues_raf_highexp(:,i));
end

delta_species = zeros(length(nodes2_minsor),1);
for i=1:length(nodes2_minsor)
    delta_species(i) = max(allValues_raf_highexp(:,nodes2_minsor(i))) ./ max(allValues_wt_raf_highexp(:,nodes2_minsor(i)));
end

[delta_Sorted, delta_index] = sort(delta_species,'descend');
sorted_delta = delta_species(delta_index);
high_delta = sorted_delta(sorted_delta > 2);
high_delta_index = find(sorted_delta > 2);
figure
bar(sorted_delta);
hold on
xlim=get(gca,'xlim');
plot(xlim,[2 2],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Maximum concentration fold change';'(Mutant / WT)'},'FontWeight','bold');
title('10-fold increased in [RAF1]');
set(gca, 'xTick', 1:46,'xTickLabel',nodelabel_minsor(delta_index),'XTickLabelRotation',45);
bh = bar(1:numel(sorted_delta),diag(sorted_delta),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(k).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end




%% Parameter sensitivity analysis of fold change in pERK for KRAS mutant vs. WT RAS (Fig. 6B)
%specify number of parameter sets 
number = 3000;
%specify parameter sampling method
sampling = 'lograndom';
%set random seed
rng(1);
% Solve ODEs with 3000 parameter arrays for wild-type RAS
[my_other_params_minsor, my_other_params_plussor, my_other_params_yinit, my_total_change_params, my_y_max_total_allparams, my_y_ave_total, my_y_steady_total, output, min_sor_total_output, plus_sor_total_output] = evaluateparam_one(wanted_param_wt, params, timeSpan, yinit, number, sampling);
% Solve ODEs with 3000 parameter arrays for KRAS mutant
[my_other_params_minsor_mutant, my_other_params_plussor_mutant, my_other_params_yinit_mutant, my_total_change_params_mutant, my_y_max_total_allparams_mutant, my_y_ave_total_mutant, my_y_steady_total_mutant, output_mutant, min_sor_total_output_mutant, plus_sor_total_output_mutant] = evaluateparam_one(wanted_param_mutant, params_mutant, timeSpan, yinit_mutant, number, sampling);

%Set change in pERK as Y in PLSR
pERK_change = abs(my_y_max_total_allparams{6} - my_y_max_total_allparams_mutant{6}); % change caused by parameter variations

ncomp = 10;
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(number, sampling, my_other_params_plussor, pERK_change, ncomp, wanted_param_wt);

%plot_plsr(PLSR_cat_mutant, number, sampling, Xloadings, vipScores, PCTVAR, 'Xloadings');
plot_plsr(PLSR_cat_mutant, number, sampling, Xloadings, vipScores, PCTVAR, 'vip');
%%
save rasmutant_main.mat
