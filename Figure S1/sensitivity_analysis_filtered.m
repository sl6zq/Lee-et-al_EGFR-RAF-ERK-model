clear all
clc
% setting to save ODE solutions for dependent variable in PLSR SA 
savedata = false;
% This script plots 
% B) Comparison of PLSR PC loadings onto the first 6 PCs after removing parameter sets that produce large aggregate error against 6 sets of experimental data and parameter
% eigenvalues from model Sloppiness Analysis onto the 0~5 modes 
% C) Average difference in parameter loadings compared to previous set
% D) Q2Y of PLSR SA with filtered parameter sets
% E) Jaccard index plots that compares similarity among unique parameters
% F) Pearson's correlation coefficients between scaled sensitivity values from PLSR SA and Sloppiness analysis 
% G) Comparison of random vs. LHS in PLSR SA

%% Define plot color scheme and save settings
cscheme_vectors = {[90/255 0/255 160/255; 1 1 1; 235/255, 107/255, 35/255]; [148/255 103/255 186/255; 1 1 1; 235/255, 107/255, 35/255]; [113/255 35/255 205/255; 1 1 1; 235/255, 107/255, 35/255]; [124/255 84/255 162/255; 1 1 1; 235/255, 107/255, 35/255]; [25 23 140; 255 255 255; 201 69 122]./255};
final_cscheme   = cscheme_vectors{end};
% save settings
if savedata == 1
    disp('Data and plots will be saved')
    CurrentDir = pwd;
    strnow                  = [datestr(now, 'dd-mmm-yyyy-HH-MM')];
    % Check MATLAB/Windows version compatibility
    if ispc
        delim_ch = '\';
    else
        delim_ch = '/';
    end
    mkdir([CurrentDir delim_ch 'Plots_generated' delim_ch strnow delim_ch]);
    folder = [CurrentDir delim_ch 'Plots_generated' delim_ch strnow delim_ch];
else
    disp('Data and plots will NOT be saved')
end
%% Define parameters
rng(0);    %random seed for SA
paramlist; %baseline parameters (best fit param values from SloppyCell)
sampling         = 'random'; %sampling method for SA, choices: random, LHS, log-random, log-LHS
% ODE param
num_Param        = 3000; % at least 2 for cv partition, num_Param >= 122
timeSpan         = 0:1:60;

% Retain parameter sets that have cost values in top 1/3
remove_threshold = round(num_Param * 1/3);

% PLSR model with 40 componenets (max. possible # of components)
ncomp            = 40;

%% B: Comparison of PLSR SA with filtered parameter sets vs. Sloppiness analysis
% Get ODE solutions with 3000 parameter sets 
% Solve ODEs
lowraf1_out_porder = evaluateparam_one(wanted_param_porder, params, timeSpan, yinit, num_Param, sampling);
if savedata == 1
    save ODE_slns_3000.mat
end
% Filter parameter sets
% Remove parameters that are not within error
if (remove_threshold < ncomp)
%    remove_threshold = ncomp+1;
   disp("1/3 of filtered param set is less than ncomp");   
   return
end

[my_remain_params, my_y_max_remain, aggregate_error] = filter_fnc(num_Param, remove_threshold, timeSpan, lowraf1_out_porder.my_other_params_minsor, lowraf1_out_porder.min_sor_total_output, lowraf1_out_porder.plus_sor_total_output, lowraf1_out_porder.my_y_max_total, 'off');
datasets;
% Fit model with filtered params
agg_err_pred      = zeros(num_Param,1);
min_sor_ras_fits  = zeros(num_Param,length(min_sor_rastimedata));
plus_sor_ras_fits = zeros(num_Param,length(plus_sor_rastimedata));
min_sor_raf_fits  = zeros(num_Param,length(min_sor_raftimedata));
plus_sor_raf_fits = zeros(num_Param,length(plus_sor_raftimedata));
min_sor_perk_fits = zeros(num_Param,length(min_sor_perktimedata));
min_sor_pmek_fits = zeros(num_Param,length(min_sor_pmektimedata));
%calculate cost for each param set
for i =1:num_Param
    min_sor_ras_fits(i,:)  = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),min_sor_rastimedata,'min_sor_rastimedata', wanted_param_porder);
    plus_sor_ras_fits(i,:) = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),plus_sor_rastimedata,'plus_sor_rastimedata', wanted_param_porder);
    min_sor_raf_fits(i,:)  = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),min_sor_raftimedata,'min_sor_raftimedata', wanted_param_porder);
    plus_sor_raf_fits(i,:) = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),plus_sor_raftimedata,'plus_sor_raftimedata', wanted_param_porder);
    min_sor_perk_fits(i,:) = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),min_sor_perktimedata,'min_sor_perktimedata', wanted_param_porder);
    min_sor_pmek_fits(i,:) = fullEGFR9_onemodel_fit_40_v2(lowraf1_out_porder.my_other_params_minsor(i,:),min_sor_pmektimedata,'min_sor_pmektimedata', wanted_param_porder);
end

for j=1:3000
    agg_err_pred(j) = 1/2 * (sum((min_sor_ras_fits(j,:)' - min_sor_rasdatanums).^2./(min_sor_rasdataerror.^2))     ...
                           + sum((min_sor_raf_fits(j,:)' - min_sor_rafdatanums).^2./(min_sor_rafdataerror.^2))  ...
                           + sum((min_sor_perk_fits(j,:)' - min_sor_perkdatanums).^2./(min_sor_perkdataerror.^2)) ...
                           + sum((min_sor_pmek_fits(j,:)' - min_sor_pmekdatanums).^2./(min_sor_pmekdataerror.^2)) ...
                           + sum((plus_sor_ras_fits(j,:)' - plus_sor_rasdatanums).^2./(plus_sor_rasdataerror.^2)) ...
                           + sum((plus_sor_raf_fits(j,:)' - plus_sor_rafdatanums).^2./(plus_sor_rafdataerror.^2)));
end


% PLSR with updated param set
[plsr_z_score_x_remain, plsr_z_score_y_remain, Xloadings_remain, vipScores_remain, Q2Y_remain, R2Y_remain, BETA_remain, PCTVAR_remain, PC1loadings_remain] = plsr_fnc(remove_threshold, sampling, my_remain_params, my_y_max_remain, ncomp, wanted_param);
names2 = {'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
names_vec2 = ["kEf" "kcatE" "kpMEK" "knfpSOS" "kBr" "kBf" "kdpERK" "kR1r" "kR1f" "knfpBR1" "kpR1" "kSOSr" "kSOSf" "kRgneslow" "kSon" "kG2SOSr" "kG2SOSf" "kiR1r" "kiR1f" "kSoff" "kpERK" "kRhydro" "Kmgneslow" "kdpMEK" "knfpiR1r" "kdR1r" "kiBr" "kdR1f" "kiBf" "kfpBr" "kdpSOS" "kdp" "kdEr" "kdE,f" "kdpR1" "kG2r" "kG2f" "knfpiBr" "kfpR1r" "kEr"]';
PLSR_cat2 = categorical(names2);
PLSR_cat2 = reordercats(PLSR_cat2,{'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'});
% match the pathway order
wanted_param_porder = [2;163;129;127;5;120;147;143;53;50;67;66;10;114;...
57;102;95;27;25;107;105;162;49;130;22;16;45;111;159;65;77;75;71;108;106;...
104;8;103;86;23];
% Make Xloadings follow pathway order
[Xloadings_remain_ordered] = order_data(Xloadings_remain, PLSR_cat2, PLSR_cat_porder);
[norm_loadings_remain, data_matlab_remain, data_sloppy_remain] = sloppy_comparison(remove_threshold, sampling, ncomp, select_ncomp, wanted_param, Xloadings_remain_ordered, Q2Y_remain, PLSR_cat_porder, final_cscheme);

if savedata == 1
    strnow                  = [datestr(now, 'dd-mmm-yyyy-HH-MM')];
    filename = (['PLSR_Sloppy_correlation_filteredparams_' strnow '.pdf']);
    saveas(gcf, filename);
end

[threshold, jaccard_remainparams] = plot_jaccard(data_matlab_remain, data_sloppy_remain, names, names_vec);
if savedata == 1
    saveas(gcf,[pwd '/Plots/PLSR_Sloppy_jaccardidx_filteredparams.pdf']);
end

%% C: Average difference in parameter loadings compared to previous parameter set
paramlist;
timeSpan = 0:1:60;
ncomp = 40;

n1 = 50;
n2 = 100;
n3 = 3000;
x = n1:n2:n3;

listplsr_random      = zeros(n3/n2,1);
listq2y_random       = zeros(n3/n2,1);
iteration            = 0;
plsr_sens_old_random = zeros(40,1);
listimportant_random = zeros(n3/n2,1);
for num = n1:n2:n3
    iteration = iteration + 1
    output    = evaluateparam_one(wanted_param, params, timeSpan, yinit, num, 'random');
    [plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, PCTVAR, PC1loadings] = plsr_fnc(num, 'random', output.my_other_params_plussor, output.my_y_max_total_allparams, ncomp, wanted_param);
    
    PC1loadings_random = PC1loadings(:,1);
    Q2Y_random = Q2Y;
    
    changeplsr_random = abs(PC1loadings_random - plsr_sens_old_random);
    plsr_sens_old_random = PC1loadings_random;

    ave_change_plsr_random = mean(changeplsr_random);
   
    listplsr_random(iteration) = ave_change_plsr_random;
    listq2y_random(iteration,:) = Q2Y_random(:,30);
   
    clear ave_change_plsr_random
    clear Q2Y_random
end

fig = figure; 
plot(x,listplsr_random,'or','MarkerFaceColor','r')
fig.Position = [100 100 450 330]
xlabel('Number of random parameter sets', 'FontWeight','bold')
ylabel({'Mean change in PLSR parameter loadings'; 'compared to previous parameter set'}, 'FontWeight','bold')
set(gcf,'color','w')
set(gca,'FontSize',8)
hold on
hold off

%% D: Q2Y of PLSR SA with filtered parameter sets
figure;
plot(1:ncomp,100*Q2Y_remain,'-bo');
xlabel('Components')
ylabel('Q^2Y')
ylim([0 100])
%% E: Plot comparison between original PLSR jaccard idx vs updated PLSR jaccard idx 
labels_PCs = {['PC 1'], ['PC 1~2'], ['PC 1~3'], ['PC 1~4'], ['PC 1~5'], ['PC 1~6']};
h = figure;
for i=1:6
    subplot(3,2,i);
    plot(threshold, jaccard_remainparams(:,i),'-o','MarkerFaceColor', 'b');
    hold on
    xlabel('Threshold','FontSize',8);
    ylabel('Jaccard index','FontSize',8);
    plot(threshold, jaccard_allparams(:,i),'-o','MarkerFaceColor', 'r');
    legend('1000 filtered parameter sets','3000 parameter sets');
    title(labels_PCs{i}, 'FontSize',8);
    filename = ['jaccardcomparison' num2str(i) '.fig'];
    %saveas(h,filename);
end
hold off

labels_PCs = {['PC 1'], ['PC 1~2'], ['PC 1~3'], ['PC 1~4'], ['PC 1~5'], ['PC 1~6']};
for i=1:2
    h = figure;
    h.Position = [100 100 450 350];
    plot(threshold, jaccard_remainparams(:,i),'-x','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','LineWidth',0.5);
    hold on
    xlabel('Threshold','FontSize',10, 'FontName', 'Arial');
    ylabel('Jaccard index','FontSize',10, 'FontName', 'Arial');
    plot(threshold, jaccard_allparams(:,i),'-o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r','LineWidth',0.5);
    legend('1000 filtered parameter sets','3000 parameter sets', 'FontSize',10, 'FontName', 'Arial');
    title(labels_PCs{i}, 'FontSize',10, 'FontName', 'Arial');
end
hold off
%% F: Pearson's correlation coefficients between scaled sensitivity values from PLSR SA and Sloppiness analysis 
plot_correlation(remove_threshold, sampling, abs(data_matlab_remain), abs(data_sloppy_remain), names);
%% G: Compare sampling methods 
sampling2   = 'LHS'; 
lowraf1_out = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);
lowraf1_out_LHS = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling2);
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(num_Param, sampling, lowraf1_out.my_other_params_plussor, lowraf1_out.my_y_max_total, ncomp, wanted_param);
[plsr_z_score_x_LHS, plsr_z_score_y_LHS, Xloadings_LHS, vipScores_LHS, Q2Y_LHS, R2Y_LHS, BETA_LHS, PCTVAR_LHS, PC1loadings_LHS] = plsr_fnc(num_Param, sampling2, lowraf1_out_LHS.my_other_params_plussor, lowraf1_out_LHS.my_y_max_total, ncomp, wanted_param);
figure
loadings_comb = [Xloadings(:,1) Xloadings_LHS(:,1)];
bar(PLSR_cat, loadings_comb);
t0 = {['PLSR Principal Component 1 (' num2str(num_Param) ' ' sampling ' parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
title(t0);
ylabel('Parameter loadings','FontWeight','bold','FontSize',8);
set(gcf,'color','w');
legend('Random', 'LHS');
ax = gca; ax.TitleFontSizeMultiplier = 1;
set(gca,'FontSize',8);
if savedata == 1
    strnow   = [datestr(now, 'dd-mmm-yyyy-HH-MM')];
    filename = (['random_LHS_PLSRSA_' strnow]);
    file     = fullfile(folder,filename);
    saveas(gcf, file, 'pdf')
end
