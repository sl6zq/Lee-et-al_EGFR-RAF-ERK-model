clear all
clc
% If savedata is true, save all plots and data in .mat files
savedata = false;
% 1) 4A: Comparison of PLSR PC loadings onto the first 6 PCs and parameter Parameter sensitivity plots using PLSR - loadings
% 2) 4B: Pearson's correlation coefficients between rescaled parameter loadings (PLSR SA) and Sloppiness analysis eigenvalues
% 3) 4C: Eigenvalues from model Sloppiness Analysis onto the 0~5 modes 
% 4) 4E: Parameters that control EGFR-ERK signaling dynamics for different RAF1 abundances
% 5) 4F: PLSR coefficents and univariate sensitivities 

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

%% Load best-fit parameters
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

%% Get ODE solutions with 3000 parameter sets 
% Solve ODEs
lowraf1_out         = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);

%Sanity check: outputs from evaluateparam_one is consistent?
check_max = [lowraf1_out.my_y_max_total{1,1}(1), lowraf1_out.my_y_max_total{2,1}(1), lowraf1_out.my_y_max_total{3,1}(1), lowraf1_out.my_y_max_total{4,1}(1), lowraf1_out.my_y_max_total{5,1}(1), lowraf1_out.my_y_max_total{6,1}(1)];
disp(check_max)
check_t_int = [lowraf1_out.t_int_total{1,1}(1), lowraf1_out.t_int_total{2,1}(1), lowraf1_out.t_int_total{3,1}(1), lowraf1_out.t_int_total{4,1}(1), lowraf1_out.t_int_total{5,1}(1), lowraf1_out.t_int_total{6,1}(1)];
disp(check_t_int)

if savedata == 1
    save ODE_slns_3000.mat
end
%% 4A: PLSR model (low RAF1 abundance)
% Select number of components to compare with Sloppy analysis (e.g., 1-6)
select_ncomp       = {'1','2','3', '4', '5', '6'};
lowraf1_out_porder = evaluateparam_one(wanted_param_porder, params, timeSpan, yinit, num_Param, sampling);
% Comparing loadings between PLSR and Sloppy 
[plsr_z_score_x, plsr_z_score_y, Xloadings_porder, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(num_Param, sampling, lowraf1_out_porder.my_other_params_plussor, lowraf1_out_porder.my_y_max_total, ncomp, wanted_param_porder);
names2 = {'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
names_vec2 = ["kEf" "kcatE" "kpMEK" "knfpSOS" "kBr" "kBf" "kdpERK" "kR1r" "kR1f" "knfpBR1" "kpR1" "kSOSr" "kSOSf" "kRgneslow" "kSon" "kG2SOSr" "kG2SOSf" "kiR1r" "kiR1f" "kSoff" "kpERK" "kRhydro" "Kmgneslow" "kdpMEK" "knfpiR1r" "kdR1r" "kiBr" "kdR1f" "kiBf" "kfpBr" "kdpSOS" "kdp" "kdEr" "kdE,f" "kdpR1" "kG2r" "kG2f" "knfpiBr" "kfpR1r" "kEr"]';
PLSR_cat2  = categorical(names2);
PLSR_cat2  = reordercats(PLSR_cat2,{'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'});
[Xloadings_porder] = order_data(Xloadings, PLSR_cat2, PLSR_cat_porder);
[norm_loadings, data_matlab, data_sloppy] = sloppy_comparison(num_Param, sampling, ncomp, select_ncomp, wanted_param_porder, Xloadings_porder, Q2Y, PLSR_cat_porder, final_cscheme);
%% 4B: Percent variance explained by each mode/PLS Component in PLSR SA & Sloppiness analysis
%load eigenvalues of each mode, obtained from SloppyCell 
eigvaldata = load('egfr-logsingval.mat');
raw_vals   = zeros(length(eigvaldata.logs),1);
cumulative = 0;
cumulativesums = zeros(length(eigvaldata.logs),1);
for i=1:length(eigvaldata.logs)
    raw_vals(i) = exp(-i);
    cumulative =  cumulative + raw_vals(i);
    cumulativesums(i) = cumulative;
end
all_modes_eigsum = sum(cumulative);
percent_eigval   = cumulativesums./all_modes_eigsum * 100;

inf_1to6 = zeros(6,1);
for i=1:6 
    inf_1to6(i) = exp(eigvaldata.logs(i));
end
inf_1to6_sum = sum(inf_1to6);

inf_1to2 = zeros(2,1);
for i=1:2 
    inf_1to2(i) = exp(eigvaldata.logs(i));
end
inf_1to2_sum = sum(inf_1to2);

inf_1to41 = zeros(41,1);
for k=1:40 
    inf_1to41(k) = exp(eigvaldata.logs(k));
end
inf_1to41_sum = sum(inf_1to41);

inf_6modes = inf_1to6_sum/inf_1to41_sum * 100;
stment = ['6 stiffest modes explain ' num2str(inf_6modes) '% of model behavior'];
disp(stment)
inf_2modes = inf_1to2_sum/inf_1to41_sum * 100;
stment2 = ['2 stiffest modes explain ' num2str(inf_2modes) '% of model behavior'];
disp(stment2)

figsize = [100 30 300 200];
h6 = figure; h6.Position = figsize;
yyaxis left
plot(1:40, eigvaldata.logs,'ko','MarkerSize',4);
ylabel('log(\lambda)');
ylim([min(eigvaldata.logs)*1.1,max(eigvaldata.logs)*1.5]);
yyaxis right
plot(1:40, percent_eigval,'k*','MarkerSize',4);
ylim([0,110]); ylabel('% fit to data explained');
xlabel('Mode');
%% 4C: Correlation coefficients between loadings from PLSR with sloppy analysis eigenvalues onto each mode (Fig. 3A)
plot_correlation(num_Param, sampling, abs(data_matlab), abs(data_sloppy), names_porder); %names are in same order for PLSR SA & Sloppy
if savedata == 1
    saveas(gcf,[pwd '/Plots/PLSR_Sloppy_correlation.pdf']);
end
%% 4E: Parameters that control EGFR-ERK signaling dynamics for different RAF1 abundances
mult = [0.01 0.1 1 10 100];
for i=1:length(mult)
    yinit_highraf1 = yinit;
    yinit_highraf1(15) = yinit_highraf1(15) * mult(i);
    title_input = ['[RAF1] * ' num2str(mult(i))];
    diffraf1_out   = evaluateparam_one(wanted_param, params, timeSpan, yinit_highraf1, num_Param, sampling);
    % PLSR w/ maximum outputs: +/-sor RASGTP, +/-memb RAF1, pMEK, pERK
    [diffraf1_plsr_z_score_x, diffraf1_plsr_z_score_y, diffraf1_Xloadings, diffraf1_vipScores, diffraf1_Q2Y, diffraf1_R2Y, diffraf1_BETA, diffraf1_PCTVAR, diffraf1_PC1loadings] = plsr_fnc(num_Param, sampling, diffraf1_out.my_other_params_plussor, diffraf1_out.my_y_max_total, ncomp, wanted_param);
    plot_plsr(PLSR_cat, num_Param, sampling, diffraf1_Xloadings, diffraf1_vipScores, diffraf1_PCTVAR, 'vip', title_input);
end
if savedata == 1
    save('diff_RAF1_vip_1xRAF1.mat', 'diffraf1_out', 'diffraf1_plsr_z_score_x', 'diffraf1_plsr_z_score_y', 'diffraf1_Xloadings', 'diffraf1_vipScores', 'diffraf1_Q2Y', 'diffraf1_R2Y', 'diffraf1_BETA', 'diffraf1_PCTVAR', 'diffraf1_PC1loadings')
end
%% 4F: Local, Univariate, PLSR SA comparison 
% Univariate SA
load('univsens_int.mat')
% Load local SA results from Julia
local_SA                = readtable('norm_integrated_lsa_res_-sorafenib.csv');
local_SA_paramvals      = table2array(local_SA(:,1:4));
local_SA_paramvals_norm = zeros(size(local_SA_paramvals,1), size(local_SA_paramvals,2));
for k=1:size(local_SA_paramvals,2)-1
   max_col = max(local_SA_paramvals(:,k));
   local_SA_paramvals_norm(:,k) = local_SA_paramvals(:,k)./max_col;
end
local_SA_paramvals_norm(:,4) = local_SA_paramvals(:,4);
localSA_paramorder     = table2cell(local_SA(:,5));
% Edit kpMEK so that it is in same text format as other parameters
localSA_paramorder(29) = {'kpMEK'};
wanted_localSA_paramorder = [2;163;129;127;5;120;147;143;67;66;53;...
50;57;102;95;27;25;162;22;16;45;107;105;49;130;103;23;114;8;86;45;10];
lowraf1_out_localSAporder = evaluateparam_one(wanted_localSA_paramorder, params, timeSpan, yinit, num_Param, sampling);
% Comparing loadings between PLSR and Sloppy 
ncomp_localSA = length(wanted_localSA_paramorder);
[plsr_z_score_x_localSAporder, plsr_z_score_y_localSAporder , Xloadings_localSAporder, vipScores_localSAporder , Q2Y_localSAporder , R2Y_localSAporder , BETA_localSAporder , PCTVAR_localSAporder, PC1loadings_localSAporder ] = plsr_fnc(num_Param, sampling, lowraf1_out_localSAporder.my_other_params_plussor, lowraf1_out_localSAporder.my_y_max_total, ncomp_localSA, wanted_localSA_paramorder);

coeff_percentmax_int_porder = zeros(length(wanted_localSA_paramorder), 1);
for i=1:6
    coeff_percentmax_int_porder(:,i) = abs(BETA_localSAporder(2:end,i)');
end
coeff_percentmax_int_comb_porder = sum(coeff_percentmax_int_porder,2);
coeff_percentmax_int_comb_porder = coeff_percentmax_int_comb_porder/max(coeff_percentmax_int_comb_porder) * 100;

local_SA_paramvals_comb = sum(local_SA_paramvals,2);
local_SA_paramvals_comb = local_SA_paramvals_comb/max(local_SA_paramvals_comb) * 100;

% Concatenate 3 matrices (dim, A1, A2, A3)
sens_all_int    = cat(2,local_SA_paramvals_comb, all_local_outputs', coeff_percentmax_int_comb_porder);
all_species_int = {'Local', 'Univariate', 'Multivariate'};
% Generate plot
f = figure;
f.Position = [100 100 800 130];
h1 = heatmap(localSA_paramorder', all_species_int, sens_all_int','CellLabelColor','none', 'GridVisible','off','FontSize',9);
ax = gca;
set(gca,'FontSize',9);
bwr = @(n)interp1([1 2 3], final_cscheme, linspace(1, 3, n), 'linear');
colormap(bwr(200)); 
if savedata == 1
    strnow   = [datestr(now, 'dd-mmm-yyyy-HH-MM')];
    filename = (['uni_local_PLSRSA_comparison_' strnow]);
    file     = fullfile(folder,filename);
    saveas(gcf, file, 'pdf')
end