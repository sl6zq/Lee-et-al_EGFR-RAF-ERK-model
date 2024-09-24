clear all
clc
% If savedata is true, save all plots and data in .mat files
savedata = false;
% 1) A: cost, computational time over iterations for all vs. VIP > 1 params varied
% 2) B: RAS-GTP, membrane RAF1, pMEK, pERK dynamics when all vs. VIP > 1 params are varied 
% 3) C: cost and wall time comparison between VIP > 1 vs. all parameters
% varied
% 3) D: pMEK and pERK compared between all vs. VIP > 1 params
% 4) E: sensitivity of RAS-GTP to sloppy vs. stiff params varied
% 5) F: cost comparison between model fitting by varying all vs. subset of
% params
% 6) G: sensitivity of model fit to individual data point from Sloppiness analysis

%% Load exp data
% Multiple optimization runs (number of runs = datasample)
datasample = 20;
% Load experimental data (6 sets)
datasets;

%% Define parameters
rng(1) %random seed for SA
paramlist; % load best-fit values from SloppyCell
%baselineparam_init_list; %baseline parameters 
sampling  = 'random'; %sampling method for SA

% ODE param
num_Param = 3000; % at least 2 for cv partition, num_Param >= 122
timeSpan  = 0:1:60;
% Maximum number of PLS components in PLSR model
ncomp     = 40;

%% A
% cost, computational time over iterations for all vs. VIP > 1 params varied
baselineparam_init_list; % Initial point for model fitting: baseline parameter values 
all_my_wanted_params = params(wanted_param)'; % Only get wanted parameters (model rate constants)
maxfuncevals         = 40000;
maxiter              = 40000;
objfunc_single_all = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata',wanted_param, yinit, params) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata',wanted_param, yinit, params) - min_sor_rafdatanums).^2./ (min_sor_rafdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param, yinit, params) - plus_sor_rasdatanums).^2./ (plus_sor_rasdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param, yinit, params) - plus_sor_rafdatanums).^2./ (plus_sor_rafdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata',wanted_param, yinit, params) - min_sor_perkdatanums).^2./ (min_sor_perkdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param, yinit, params) - min_sor_pmekdatanums).^2./ (min_sor_pmekdataerror).^2));

LB_single_all = zeros(1,length(wanted_param))';
tic
[fittedparams_single_all, fval_singleall, fval_singleallhistory, obj_val_his_singleall] = hist_fminsearch(all_my_wanted_params, LB_single_all, objfunc_single_all, maxfuncevals, maxiter);
elaspsedtime_all = toc;

% Get a subset of parameters to vary via PLSR SA
lowraf1_out = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);% Solve ODEs
% Build PLSR 
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, PCTVAR, PC1loadings] = plsr_fnc(num_Param, sampling, lowraf1_out.my_other_params_plussor, lowraf1_out.my_y_max_total, ncomp, wanted_param);
% Find parameters w/ VIP > 1
vip_wanted_param_values  = find(vipScores>1);
vip_wanted_param_index   = wanted_param(vip_wanted_param_values);
wanted_param2            = vip_wanted_param_index;
my_wanted_params         = params(wanted_param2)';
LB_single                = zeros(1,length(wanted_param2))';

objfunc_single = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata',wanted_param2, yinit, params) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata',wanted_param2, yinit, params) - min_sor_rafdatanums).^2./ (min_sor_rafdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param2, yinit, params) - plus_sor_rasdatanums).^2./ (plus_sor_rasdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param2, yinit, params) - plus_sor_rafdatanums).^2./ (plus_sor_rafdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata',wanted_param2, yinit, params) - min_sor_perkdatanums).^2./ (min_sor_perkdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param2, yinit, params) - min_sor_pmekdatanums).^2./ (min_sor_pmekdataerror).^2));

tic
[fittedparams_single, fval_singlevip, fval_singleviphistory, obj_val_his] = hist_fminsearch(my_wanted_params, LB_single, objfunc_single, maxfuncevals, maxiter);
elaspsedtime_subset = toc;

f1 = figure;
f1.Position = [100 100 400 250];
iter_subset = 0:1:length(obj_val_his)-1;
iter_all    = 0:1:length(obj_val_his_singleall)-1;
plot(iter_subset, obj_val_his, 'o', 'MarkerFaceColor', [253/255, 136/255, 45/255], 'MarkerEdgeColor',[253/255, 136/255, 45/255], 'MarkerSize',4)
hold on
plot(iter_all, obj_val_his_singleall,'o', 'MarkerFaceColor', [43/255, 133/255, 39/255], 'MarkerEdgeColor',[43/255, 133/255, 39/255], 'MarkerSize', 4);
xline(iter_subset(end), '--k', '1.33 hr');
xline(iter_all(end), '--k', '4.74 hr');
ylab = ['C (' char(952) ')'];
xlabel('Iteration'); ylabel(ylab);
legend('VIP > 1', 'All','','Location','northeastoutside')
ylim([430 max(obj_val_his)+50]); xlim([0 iter_all(end)*1.2]);

%% B
% Single best fits using all vs. VIP > 1 parameters
[min_sor_rasoutput_single,T, min_sor_ras_allValues_single, param_all_single, param_varied_single] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'min_sor_rastimedata',wanted_param, yinit, params);
[plus_sor_rasoutput_single,~, plus_sor_ras_allValues_single, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'plus_sor_rastimedata',wanted_param, yinit, params);
[min_sor_rafoutput_single,~, min_sor_raf_allValues_single, ~, ~]                                  = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'min_sor_raftimedata',wanted_param, yinit, params);
[plus_sor_rafoutput_single,~, plus_sor_raf_allValues_single, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'plus_sor_raftimedata',wanted_param, yinit, params);
[min_sor_perkoutput_single,~, min_sor_perk_allValues_single, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'min_sor_perktimedata',wanted_param, yinit, params);
[min_sor_pmekoutput_single,~, min_sor_pmek_allValues_single, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_all, linspace(0,60),'min_sor_pmektimedata',wanted_param, yinit, params);

single_allparams_outputs{1} = min_sor_rasoutput_single;
single_allparams_outputs{2} = plus_sor_rasoutput_single;
single_allparams_outputs{3} = min_sor_rafoutput_single;
single_allparams_outputs{4} = plus_sor_rafoutput_single;
single_allparams_outputs{5} = min_sor_perkoutput_single;
single_allparams_outputs{6} = min_sor_pmekoutput_single;


allparams_highVIPfitted = params;
for i=1:length(wanted_param2)
    allparams_highVIPfitted(wanted_param2(i)) = fittedparams_single(i);
end
wantedparams_highVIPfitted = allparams_highVIPfitted(wanted_param);
%allparams_highVIPfitted(wanted_param) = fittedparams_single;
[min_sor_rasoutput_single_highVIP,T, min_sor_ras_allValues_single_highVIP, param_all_single_highVIP, param_varied_single] = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'min_sor_rastimedata',wanted_param, yinit, params);
[plus_sor_rasoutput_single_highVIP,~, plus_sor_ras_allValues_single_highVIP, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'plus_sor_rastimedata',wanted_param, yinit, params);
[min_sor_rafoutput_single_highVIP,~, min_sor_raf_allValues_single_highVIP, ~, ~]                                  = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'min_sor_raftimedata',wanted_param, yinit, params);
[plus_sor_rafoutput_single_highVIP,~, plus_sor_raf_allValues_single_highVIP, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'plus_sor_raftimedata',wanted_param, yinit, params);
[min_sor_perkoutput_single_highVIP,~, min_sor_perk_allValues_single_highVIP, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'min_sor_perktimedata',wanted_param, yinit, params);
[min_sor_pmekoutput_single_highVIP,~, min_sor_pmek_allValues_single_highVIP, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(wantedparams_highVIPfitted, linspace(0,60),'min_sor_pmektimedata',wanted_param, yinit, params);

single_highVIP_outputs{1} = min_sor_rasoutput_single_highVIP;
single_highVIP_outputs{2} = plus_sor_rasoutput_single_highVIP;
single_highVIP_outputs{3} = min_sor_rafoutput_single_highVIP;
single_highVIP_outputs{4} = plus_sor_rafoutput_single_highVIP;
single_highVIP_outputs{5} = min_sor_perkoutput_single_highVIP;
single_highVIP_outputs{6} = min_sor_pmekoutput_single_highVIP;


plot_fit(single_allparams_outputs, 'no');
plot_fit(single_highVIP_outputs, 'no');
%% C,D
% Multiple good fits using all vs. subset of parameters
CI_wanted  = '1SD'; % Specify wanted CI level
options    = optimset('MaxFunEvals', maxfuncevals, 'MaxIter', maxiter, 'PlotFcns', @optimplotfval);

maxiterations = cell(length(datasample));

iteration = 0;
tic
parfor f = 1:datasample
    iteration = iteration+1
    
    %save objective function for each run 
    all_ci_objfunc{f} = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata', wanted_param) - min_sor_rasdata_rangenums(f,:)').^2 ./ (min_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata', wanted_param) - min_sor_rafdata_rangenums(f,:)').^2 ./ (min_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata', wanted_param) - plus_sor_rasdata_rangenums(f,:)').^2 ./ (plus_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata', wanted_param) - plus_sor_rafdata_rangenums(f,:)').^2 ./ (plus_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata', wanted_param) - min_sor_perkdata_rangenums(f,:)').^2 ./ (min_sor_perkdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata', wanted_param) - min_sor_pmekdata_rangenums(f,:)').^2 ./ (min_sor_pmekdataerror).^2));
    
    LB = zeros(1,length(wanted_param))';
    [all_ci_fittedparams(f,:),all_ci_objfuncval(f,:), ~, all_output] = fminsearchbnd(all_ci_objfunc{f},all_my_wanted_params,LB,[],options);
    maxiterations{f} = all_output.iterations;
   
    [all_ci_min_sor_rasoutput(:,f),T(:,f), ~, all_ci_param1(:,f), all_ci_param_vary1(:,f)] = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'min_sor_rastimedata', wanted_param);
    [all_ci_plus_sor_rasoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'plus_sor_rastimedata', wanted_param);
    [all_ci_min_sor_rafoutput(:,f),~, ~, ~]  = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'min_sor_raftimedata', wanted_param);
    [all_ci_plus_sor_rafoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'plus_sor_raftimedata', wanted_param);
    [all_ci_min_sor_pmekoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'min_sor_pmektimedata', wanted_param);
    [all_ci_min_sor_perkoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(f,:), linspace(0,60),'min_sor_perktimedata', wanted_param);
    
end
all_ci_elapsedtime = toc;
%Get minimum cost and its index
[all_ci_min_objfuncval,all_ci_min_index] = min(all_ci_objfuncval);
mean_ci_cost                             = mean(all_ci_objfuncval);
sd_ci_cost                               = std(all_ci_objfuncval);
mean_disp = ['Cost by fitting all parameters = ' num2str(mean_ci_cost) ' ' char(177) ' ' num2str(sd_ci_cost)];
disp(mean_disp)

all_ci_totaloutput    = cell(6,1);
all_ci_totaloutput{1} = all_ci_min_sor_rasoutput;
all_ci_totaloutput{2} = all_ci_plus_sor_rasoutput;
all_ci_totaloutput{3} = all_ci_min_sor_rafoutput;
all_ci_totaloutput{4} = all_ci_plus_sor_rafoutput;
all_ci_totaloutput{5} = all_ci_min_sor_perkoutput;
all_ci_totaloutput{6} = all_ci_min_sor_pmekoutput;
% find best-fit (min. cost)
all_ci_min_sor_ras_optimalfit  = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'min_sor_rastimedata', wanted_param);
all_ci_plus_sor_ras_optimalfit = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'plus_sor_rastimedata', wanted_param);
all_ci_min_sor_raf_optimalfit  = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'min_sor_raftimedata', wanted_param);
all_ci_plus_sor_raf_optimalfit = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'plus_sor_raftimedata', wanted_param);
all_ci_min_sor_pmek_optimalfit = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'min_sor_pmektimedata', wanted_param);
all_ci_min_sor_perk_optimalfit = fullEGFR9_onemodel_fit_40_v2(all_ci_fittedparams(all_ci_min_index,:),linspace(0,60),'min_sor_perktimedata', wanted_param);

all_opt_fits{1} = all_ci_min_sor_ras_optimalfit;
all_opt_fits{2} = all_ci_plus_sor_ras_optimalfit;
all_opt_fits{3} = all_ci_min_sor_raf_optimalfit;
all_opt_fits{4} = all_ci_plus_sor_raf_optimalfit;
all_opt_fits{5} = all_ci_min_sor_perk_optimalfit;
all_opt_fits{6} = all_ci_min_sor_pmek_optimalfit;


time                 = T(:,1)';
all_ci_timeconf      = [time time(end:-1:1)];
all_ci_output_tscore = cell(6,1); 
all_ci_output_CI     = cell(6,1); 
all_ci_output_conf   = cell(6,1);
all_ci_output_mean   = cell(6,1);
all_ci_output_SEM    = cell(6,1);
%get confidence interval for each output and bold the mean optimization run/optimization run with
%the smallest objective function
for o=1:6
    for k=1:size(all_ci_totaloutput{1},1)
        N2                               = size(all_ci_totaloutput{1},1);
        all_ci_output_mean{o}            = mean(all_ci_totaloutput{o},2);
        all_ci_output_SEM{o}             = std(all_ci_totaloutput{o},0,2)/sqrt(N2);
        if strcmp(CI_wanted, '1SD')
            all_ci_output_tscore{o}(k,:) = tinv([0.16 0.84], N2-1); %get tscore for 68% (1 SD) CI
        elseif strcmp(CI_wanted, '2SD')
            all_ci_output_tscore{o}(k,:) = tinv([0.025 0.975], N2-1); %get tscore for 95% (2 SD) CI
        end
        all_ci_output_CI{o}(k,:)         = all_ci_output_mean{o}(k,:) + all_ci_output_tscore{o}(k,:) * all_ci_output_SEM{o}(k,:);
        all_ci_output_conf{o}            = [all_ci_output_CI{o}(:,1)' fliplr(all_ci_output_CI{o}(:,2)')];
    end
end
% Plot mean fits with CIs
plot_fit_CI(all_ci_timeconf, all_ci_output_conf, all_ci_output_mean)
% Plot all individual runs
plot_indiv_runs(all_ci_totaloutput, all_opt_fits)

% Solve for x (list of parameters) that minimize this objective function using fminsearch 
maxitr_subset = cell(length(datasample));

itr_subset = 0;
tic
parfor f = 1:datasample
    itr_subset = itr_subset+1
    
    %save objective function for each run 
    objfunc{f} = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata', wanted_param2) - min_sor_rasdata_rangenums(f,:)').^2 ./ (min_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata', wanted_param2) - min_sor_rafdata_rangenums(f,:)').^2 ./ (min_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata', wanted_param2) - plus_sor_rasdata_rangenums(f,:)').^2 ./ (plus_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata', wanted_param2) - plus_sor_rafdata_rangenums(f,:)').^2 ./ (plus_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata', wanted_param2) - min_sor_perkdata_rangenums(f,:)').^2 ./ (min_sor_perkdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata', wanted_param2) - min_sor_pmekdata_rangenums(f,:)').^2 ./ (min_sor_pmekdataerror).^2));
    
    LB = zeros(1,length(wanted_param2))';
    [fittedparams(f,:),objfuncval(f,:), ~, output] = fminsearchbnd(objfunc{f},my_wanted_params,LB,[],options);
    maxitr_subset{f} = output.iterations;
   
    [min_sor_rasoutput(:,f),T(:,f), ~, param1(:,f), param_vary1(:,f)] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_rastimedata', wanted_param2);
    [plus_sor_rasoutput(:,f),~, ~, ~]                                 = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'plus_sor_rastimedata', wanted_param2);
    [min_sor_rafoutput(:,f),~, ~, ~]                                  = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_raftimedata', wanted_param2);
    [plus_sor_rafoutput(:,f),~, ~, ~]                                 = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'plus_sor_raftimedata', wanted_param2);
    [min_sor_pmekoutput(:,f),~, ~, ~]                                 = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_pmektimedata', wanted_param2);  
    [min_sor_perkoutput(:,f),~, ~, ~]                                 = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_perktimedata', wanted_param2);
    
end
elapsedtime = toc;
%Get minimum objective function value and its index
[min_objfuncval,min_index] = min(objfuncval);
mean_ci_subset_cost = mean(objfuncval);
sd_ci_subset_cost   = std(objfuncval);
subset_mean_disp    = ['Cost by fitting VIP > 1 parameters = ' num2str(mean_ci_subset_cost) ' ' char(177) ' ' num2str(sd_ci_subset_cost)];
disp(subset_mean_disp)
% All fits
totaloutput = cell(6,1);
totaloutput{1} = min_sor_rasoutput;
totaloutput{2} = plus_sor_rasoutput;
totaloutput{3} = min_sor_rafoutput;
totaloutput{4} = plus_sor_rafoutput;
totaloutput{5} = min_sor_perkoutput;
totaloutput{6} = min_sor_pmekoutput;
% Best fits
min_sor_ras_optimalfit  = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_rastimedata', wanted_param2);
plus_sor_ras_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'plus_sor_rastimedata', wanted_param2);
min_sor_raf_optimalfit  = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_raftimedata', wanted_param2);
plus_sor_raf_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'plus_sor_raftimedata', wanted_param2);
min_sor_perk_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_perktimedata', wanted_param2);
min_sor_pmek_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_pmektimedata', wanted_param2);

time = T(:,1)';
timeconf = [time time(end:-1:1)];

output_tscore = cell(6,1); 
output_CI     = cell(6,1); 
output_conf   = cell(6,1);
output_mean   = cell(6,1);
output_SEM    = cell(6,1);
%get confidence interval for each output and bold the mean optimization run/optimization run with
%the smallest objective function

for o=1:6
    for k=1:size(totaloutput{1},1)
        N2 = size(totaloutput{1},1);
        output_mean{o} = mean(totaloutput{o},2);
        output_SEM{o}  = std(totaloutput{o},0,2)/sqrt(N2);
        if strcmp(CI_wanted, '1SD')
            output_tscore{o}(k,:) = tinv([0.16 0.84], N2-1); %get tscore for 68% (1 SD) CI
        elseif strcmp(CI_wanted, '2SD')
            output_tscore{o}(k,:) = tinv([0.025 0.975], N2-1); %get tscore for 95% (2 SD) CI
        end
        output_CI{o}(k,:)   = output_mean{o}(k,:) + output_tscore{o}(k,:) * output_SEM{o}(k,:);
        output_conf{o}        = [output_CI{o}(:,1)' fliplr(output_CI{o}(:,2)')];
    end
end

% Plot mean fits with CIs
plot_fit_CI(timeconf, output_conf, output_mean)
% Plot all individual runs (optional)
% plot_indiv_runs(all_ci_totaloutput, all_opt_fits)

%% E
% Model fits (RAS-GTP) by varying sloppy vs. stiff parameters
% Fix kG2f, kG2r or kR1f, kR1r to known values, fit the rest of the parameters
fix_params = 'kR1fkR1r'
if strcmp(fix_params, 'kG2fkG2r')
    wanted_param_fix = wanted_param;
    wanted_param_fix(37) = [];
    wanted_param_fix(36) = [];
elseif strcmp(fix_params, 'kR1fkR1r')
    wanted_param_fix = wanted_param;
    wanted_param_fix(27) = [];
    wanted_param_fix(25) = [];
end

my_wanted_params_fix = params(wanted_param_fix)';
objfunc_single_all_fix = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata',wanted_param_fix, yinit, params) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata',wanted_param_fix, yinit, params) - min_sor_rafdatanums).^2./ (min_sor_rafdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param_fix, yinit, params) - plus_sor_rasdatanums).^2./ (plus_sor_rasdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param_fix, yinit, params) - plus_sor_rafdatanums).^2./ (plus_sor_rafdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata',wanted_param_fix, yinit, params) - min_sor_perkdatanums).^2./ (min_sor_perkdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param_fix, yinit, params) - min_sor_pmekdatanums).^2./ (min_sor_pmekdataerror).^2));

LB_single_all_fix = zeros(1,length(wanted_param_fix))';
UB_single_all_fix = ones(1,length(wanted_param_fix))' * 10^2;

tic
[fittedparams_only_fix,fval_kG2fkG2r_fix,exitflag_kG2fkG2r_fix,output_kG2fkG2r_fix]  = fminsearchbnd(objfunc_single_all_fix,my_wanted_params_fix,LB_single_all_fix,UB_single_all_fix,options);
toc

% plug in fitted values for everything except kG2f, kG2r
fittedparams_single_fix                            = my_wanted_params;
for i=1:length(wanted_param_fix)
    fittedparams_single_fix(i) = fittedparams_only_fix(i);
end
% plug in known values for kG2f, kG2r
if strcmp(fix_params, 'kG2fkG2r')
    fittedparams_single_fix(37)       = 3.8*10^-4; %kG2f
    fittedparams_single_fix(36)       = 4.6*10^2; %kG2r
elseif strcmp(fix_params, 'kRf1kR1r')
    fittedparams_single_fix(27)       = 2.5 * 10^-5; %kR1f
    fittedparams_single_fix(25)       = 8.4 * 10^0; %kR1r
end
[min_sor_rasoutput_kfix,T, ~, param_all_fix, param_varied_fix] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'min_sor_rastimedata',wanted_param, yinit, params);
[plus_sor_rasoutput_kfix,~, ~, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'plus_sor_rastimedata',wanted_param, yinit, params);
[min_sor_rafoutput_kfix,~, ~, ~, ~]                                 = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'min_sor_raftimedata',wanted_param, yinit, params);
[plus_sor_rafoutput_kfix,~, ~, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'plus_sor_raftimedata',wanted_param, yinit, params);
[min_sor_perkoutput_kfix,~, ~, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'min_sor_perktimedata',wanted_param, yinit, params);
[min_sor_pmekoutput_kfix,~, ~, ~, ~]                                = fullEGFR9_onemodel_fit_40_v2(fittedparams_single_fix, linspace(0,60),'min_sor_pmektimedata',wanted_param, yinit, params);

kfix_outputs{1} = min_sor_rasoutput_kfix;
kfix_outputs{2} = plus_sor_rasoutput_kfix;
kfix_outputs{3} = min_sor_rafoutput_kfix;
kfix_outputs{4} = plus_sor_rafoutput_kfix;
kfix_outputs{5} = min_sor_perkoutput_kfix;
kfix_outputs{6} = min_sor_pmekoutput_kfix;

plot_fit(kfix_outputs, 'save',fix_params)

%% F 
% Comparison of cost between varying all, subset, stiff, sloppy parameters
f = figure;
f.Position = [100 100 200 300];
cost_comb   = [fval, 528.773, 502.709, 706.8711];
cost_labels = categorical({'All', 'VIP>1', 'kG2f, kG2r','kR1f, kR1r'})
bar(cost_labels, cost_comb);
ylabel('Cost'); xlabel('Parameters varied');

%% G 
modesensdata = load('egfr-modesens.mat');
figure
plot(1:6, modesensdata.s_mode_0(1,1:6), '-o', 'color', [0 0.6 0.3], 'LineWidth', 1.2, 'MarkerFaceColor', [0 0.6 0.3], 'MarkerSize',5);
hold on
plot(7:12, modesensdata.s_mode_0(1,29:34), '-o', 'color',[1 0.2 0.2], 'LineWidth', 1.2, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerSize',5);
plot(13:29, modesensdata.s_mode_0(1,7:23), '-o', 'color', [0 0.5 1], 'LineWidth', 1.2, 'MarkerFaceColor', [0 0.5 1], 'MarkerSize',5);
plot(30:54, modesensdata.s_mode_0(1,35:59), '-o', 'color', [1 0.5 0], 'LineWidth', 1.2, 'MarkerFaceColor', [1 0.5 0], 'MarkerSize',5);
plot(55:59, modesensdata.s_mode_0(1,60:64), '-o', 'color', [236/255 0 140/255], 'LineWidth', 1.2, 'MarkerFaceColor', [236/255 0 140/255], 'MarkerSize',5);
plot(60:64, modesensdata.s_mode_0(1,24:28), '-o', 'color', [0.2 0.8 0.8], 'LineWidth', 1.2, 'MarkerFaceColor', [0.2 0.8 0.8], 'MarkerSize',5);
plot(xlim,[0 0],'-.k','LineWidth',2);
xlabel('Data point');
ylabel('Jv');
legend('RasGTP +sor', 'RasGTP -sor', 'Raf1 +sor', 'Raf1 -sor', 'pMEK -sor', 'pERK -sor','Location', 'southeast');
if savedata == 1
    saveas(gcf,'modesens_mode0.pdf');
end

figure
plot(1:6, modesensdata.s_mode_1(1,1:6),'-o','color', [0 0.6 0.3], 'LineWidth', 1.2, 'MarkerFaceColor',[0 0.6 0.3],'MarkerSize',5);
hold on
plot(7:12, modesensdata.s_mode_1(1,29:34),'-o','color', [1 0.2 0.2], 'LineWidth', 1.2,'MarkerFaceColor',[1 0.2 0.2],'MarkerSize',5);
plot(13:29, modesensdata.s_mode_1(1,7:23),'-o','color', [0 0.5 1], 'LineWidth', 1.2,'MarkerFaceColor',[0 0.5 1],'MarkerSize',5);
plot(30:54, modesensdata.s_mode_1(1,35:59),'-o','color', [1 0.5 0], 'LineWidth', 1.2,'MarkerFaceColor',[1 0.5 0],'MarkerSize',5);
plot(55:59, modesensdata.s_mode_1(1,60:64),'-o','color', [236/255 0 140/255],'LineWidth', 1.2,'MarkerFaceColor',[236/255 0 140/255],'MarkerSize',5);
plot(60:64, modesensdata.s_mode_1(1,24:28),'-o','color', [0.2 0.8 0.8], 'LineWidth', 1.2,'MarkerFaceColor',[0.2 0.8 0.8],'MarkerSize',5);
plot(xlim,[0 0],'-.k','LineWidth',2);
xlabel('Data point');
ylabel('Jv');
ylim([-13 5]);
legend('RasGTP +sor', 'RasGTP -sor', 'Raf1 +sor', 'Raf1 -sor', 'pMEK -sor', 'pERK -sor','Location', 'southwest');
if savedata == 1
    saveas(gcf,'modesens_mode1.pdf');
end

%% 
if savedata == 1
    save modelfitting_comparison.mat
end