clear all 
clc
% This script plots 
% 1) The mean and confidence intervals (+/- 1 SD) of model fits using a
% subset of important parameters (VIP>1) from PLSR SA (not shown)
% 2) Individual optimization runs from fitting the model to synthetic data
% mapped to the experimental data (not shown)
% 3) A single best fit using a subset of important parameters (VIP>1) (Fig.
% S2, Table S5)
% 4) A single best fit using all parameters (Table S5)

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

%% Get ODE solutions with 3000 parameter sets 
% Solve ODEs
lowraf1_out = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);

%[my_other_params_minsor, my_other_params_plussor, my_other_params_yinit, my_total_change_params, my_y_max_total_allparams, my_y_ave_total, my_y_steady_total, t_int_total, output, min_sor_total_output, plus_sor_total_output] = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);
%% Use PLSR to select a subset of params to vary
% Build PLSR model
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, PCTVAR, PC1loadings] = plsr_fnc(num_Param, sampling, lowraf1_out.my_other_params_plussor, lowraf1_out.my_y_max_total, ncomp, wanted_param);
% Find parameters w/ VIP > 1
vip_wanted_param_values  = find(vipScores>1);
vip_wanted_param_index   = wanted_param(vip_wanted_param_values);
wanted_param2            = vip_wanted_param_index;
%% Load initial parameter values
baselineparam_init_list;
%% Parameter estimation with a subset of important parameters (VIP>1) that generates confidence intervals
my_wanted_params = params(wanted_param2)';
% Solve for x (list of parameters) that minimize this objective function using fminsearch 
options = optimset('MaxFunEvals', 40000, 'MaxIter', 40000, 'PlotFcns', @optimplotfval);
maxiterations = cell(length(datasample));

iteration = 0;

tic
parfor f = 1:datasample
    iteration = iteration+1
    
    %save objective function for each run 
    objfunc{f} = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata', wanted_param2) - min_sor_rasdata_rangenums(f,:)').^2 ./ (min_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata', wanted_param2) - min_sor_rafdata_rangenums(f,:)').^2 ./ (min_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata', wanted_param2) - plus_sor_rasdata_rangenums(f,:)').^2 ./ (plus_sor_rasdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata', wanted_param2) - plus_sor_rafdata_rangenums(f,:)').^2 ./ (plus_sor_rafdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata', wanted_param2) - min_sor_perkdata_rangenums(f,:)').^2 ./ (min_sor_perkdataerror).^2)...
        + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata', wanted_param2) - min_sor_pmekdata_rangenums(f,:)').^2 ./ (min_sor_pmekdataerror).^2));
    
    LB = zeros(1,length(wanted_param2))';
    [fittedparams(f,:),objfuncval(f,:), ~, output] = fminsearchbnd(objfunc{f},my_wanted_params,LB,[],options);
    maxiterations{f} = output.iterations;
   
    [min_sor_rasoutput(:,f),T(:,f), ~, param1(:,f), param_vary1(:,f)] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_rastimedata', wanted_param2);
    [plus_sor_rasoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'plus_sor_rastimedata', wanted_param2);
    [min_sor_rafoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_raftimedata', wanted_param2);
    [plus_sor_rafoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'plus_sor_raftimedata', wanted_param2);
    [min_sor_perkoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_perktimedata', wanted_param2);
    [min_sor_pmekoutput(:,f),~, ~, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams(f,:), linspace(0,60),'min_sor_pmektimedata', wanted_param2);
    
end
elapsedtime = toc
%Get minimum objective function value and its index
[min_objfuncval,min_index] = min(objfuncval);

totaloutput = cell(6,1);
totaloutput{1} = min_sor_rasoutput;
totaloutput{2} = plus_sor_rasoutput;
totaloutput{3} = min_sor_rafoutput;
totaloutput{4} = plus_sor_rafoutput;
totaloutput{5} = min_sor_perkoutput;
totaloutput{6} = min_sor_pmekoutput;

min_sor_ras_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_rastimedata', wanted_param2);
plus_sor_ras_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'plus_sor_rastimedata', wanted_param2);
min_sor_raf_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_raftimedata', wanted_param2);
plus_sor_raf_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'plus_sor_raftimedata', wanted_param2);
min_sor_perk_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_perktimedata', wanted_param2);
min_sor_pmek_optimalfit = fullEGFR9_onemodel_fit_40_v2(fittedparams(min_index,:),linspace(0,60),'min_sor_pmektimedata', wanted_param2);

time = T(:,1)';
timeconf = [time time(end:-1:1)];

output_tscore = cell(6,1); 
output_CI95 = cell(6,1); 
output_conf = cell(6,1);
output_mean = cell(6,1);
output_SEM = cell(6,1);
%get confidence interval for each output and bold the mean optimization run/optimization run with
%the smallest objective function

for o=1:6
    for k=1:size(totaloutput{1},1);
        N2 = size(totaloutput{1},1);
        output_mean{o} = mean(totaloutput{o},2);
        output_SEM{o} = std(totaloutput{o},0,2)/sqrt(N2);
        %output_tscore{o}(k,:) = tinv([0.025 0.975], N2-1); % get tscore
        %for 95% CI
        %for 95% (2 SD) CI
        output_tscore{o}(k,:) = tinv([0.16 0.84], N2-1); %get tscore for 68% (1 SD) CI
        output_CI95{o}(k,:) = output_mean{o}(k,:) + output_tscore{o}(k,:) * output_SEM{o}(k,:);
        output_conf{o} = [output_CI95{o}(:,1)' fliplr(output_CI95{o}(:,2)')];
    end
end

figure
%subplot(1,3,1);
p1 = fill(timeconf,output_conf{1}, [1 0.2 0.2],'FaceAlpha',0.3);
p1.FaceColor = [1, 0.2, 0.2];
p1.EdgeColor = 'none';
hold on
p2 = fill(timeconf,output_conf{2},[0 0.6 0.3],'FaceAlpha',0.3);
p2.FaceColor = [0 0.6 0.3];
p2.EdgeColor = 'none';
plot(min_sor_rastimedata,min_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', [1 0.2 0.2],'LineStyle', 'none');
plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', [0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
e1 = errorbar(min_sor_rastimedata,min_sor_rasdatanums,min_sor_rasdataerror,'LineStyle','none');
e1.Color = 'black';
e1.LineWidth = 1;
%p12 = plot(linspace(0,60),min_sor_ras_optimalfit,'red','LineWidth',2);
%p22 = plot(linspace(0,60), plus_sor_ras_optimalfit,'blue','LineWidth',2);
p12 = plot(linspace(0,60),output_mean{1},'Color',[1 0.2 0.2],'LineWidth',2);
p22 = plot(linspace(0,60), output_mean{2},'Color',[0 0.6 0.3],'LineWidth',2);
xlabel('Time (min)')
ylabel('Concentration (molec/cell)');
legend([p12 p22],{'-sor GTP-bound Ras' '+sor GTP-bound Ras'});

figure
p3 = fill(timeconf,output_conf{3},[1 0.5 0],'FaceAlpha',0.3);
p3.FaceColor = [1 0.5 0];
p3.EdgeColor = 'none';
hold on
p4 = fill(timeconf,output_conf{4},[0 0.5 1],'FaceAlpha',0.3);
p4.FaceColor = [0 0.5 1];
p4.EdgeColor = 'none';
plot(min_sor_raftimedata,min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0],'LineStyle', 'none');
e2 = errorbar(min_sor_raftimedata,min_sor_rafdatanums,min_sor_rafdataerror,'LineStyle','none');
e2.Color = 'black';
e2.LineWidth = 1;
plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor', [0 0.5 1],'LineStyle', 'none');
e3 = errorbar(plus_sor_raftimedata,plus_sor_rafdatanums,plus_sor_rafdataerror,'LineStyle','none');
e3.Color = 'black';
e3.LineWidth = 1;
%p32 = plot(linspace(0,60),min_sor_raf_optimalfit,'green','LineWidth',2);
%p42 = plot(linspace(0,60), plus_sor_raf_optimalfit,'magenta','LineWidth',2);
p32 = plot(linspace(0,60), output_mean{3},'Color',[1 0.5 0],'LineWidth',2);
p42 = plot(linspace(0,60), output_mean{4},'Color',[0 0.5 1],'LineWidth',2);
xlabel('Time (min)')
ylabel('Concentration (molec/cell)');
legend([p32 p42],{'-sor membrane Raf1' '+sor membrane Raf1'});
 
figure
p5 = fill(timeconf,output_conf{5},[0.2 0.8 0.8],'FaceAlpha',0.3);
p5.FaceColor = [0.2 0.8 0.8];
p5.EdgeColor = 'none';
hold on
p6 = fill(timeconf,output_conf{6},[236/255 0 140/255],'FaceAlpha',0.3);
p6.FaceColor = [236/255 0 140/255];
p6.EdgeColor = 'none';
plot(min_sor_perktimedata,min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor', [0.2 0.8 0.8], 'MarkerEdgeColor', [0.2 0.8 0.8],'LineStyle', 'none');
e4 = errorbar(min_sor_perktimedata,min_sor_perkdatanums,min_sor_perkdataerror,'LineStyle','none');
e4.Color = 'black';
e4.LineWidth = 1;
plot(min_sor_pmektimedata,min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor', [236/255 0 140/255], 'MarkerEdgeColor',[236/255 0 140/255],'LineStyle', 'none');
e5 = errorbar(min_sor_pmektimedata,min_sor_pmekdatanums,min_sor_pmekdataerror,'LineStyle','none');
e5.Color = 'black';
e5.LineWidth = 1;
%p52 = plot(linspace(0,60),min_sor_perk_optimalfit,'black','LineWidth',2);
%p62 = plot(linspace(0,60), min_sor_pmek_optimalfit,'cyan','LineWidth',2);
p52 = plot(linspace(0,60),output_mean{5},'Color',[0.2 0.8 0.8],'LineWidth',2);
p62 = plot(linspace(0,60), output_mean{6},'Color',[236/255 0 140/255],'LineWidth',2);
xlabel('Time (min)')
ylabel('Concentration (molec/cell)');
legend([p52 p62],{'-sor pERK' '-sor pMEK'});       
hold off
%% Parameter estimation with VIP>1, excluding RAF1 parameters



%% Plot individual optimization runs
f1 = figure;
f2 = figure;
f3 = figure;
for k = 1:datasample
    set(0, 'CurrentFigure', f1);
    p11 = plot(linspace(0,60),min_sor_rasoutput(:,k),'Color',[1 0.2 0.2],'LineWidth', 2);
    hold on
    p11.Color(4) = 0.25;
    p12 = plot(linspace(0,60),min_sor_ras_optimalfit,'Color',[1 0.2 0.2],'LineWidth', 2);
    p13 = plot(min_sor_rastimedata,min_sor_rasdatanums, 'Marker', 'o', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', [1 0.2 0.2], 'LineStyle', 'none');
    e1 = errorbar(min_sor_rastimedata,min_sor_rasdatanums,min_sor_rasdataerror,'LineStyle','none');
    e1.Color = 'black';
    e1.LineWidth = 1.4;
    xlabel('Time (min)');
    ylabel('Concentration (molec/cell)');
    p21 = plot(linspace(0,60),plus_sor_rasoutput(:,k),'Color',[0 0.6 0.3],'LineWidth', 2);
    p21.Color(4) = 0.25;
    p22 = plot(linspace(0,60),plus_sor_ras_optimalfit,'Color',[0 0.6 0.3],'LineWidth', 2);
    p23 = plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3],'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
    legend([p12 p22],{'-sor GTP-bound total Ras' '+sor GTP-bound total Ras'});

    set(0, 'CurrentFigure', f2);
    p31 = plot(linspace(0,60),min_sor_rafoutput(:,k),'Color',[1 0.5 0],'LineWidth', 2);
    hold on
    p31.Color(4) = 0.25;
    p32 = plot(linspace(0,60),min_sor_raf_optimalfit,'Color',[1 0.5 0],'LineWidth', 2);
    p33 = plot(min_sor_raftimedata,min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor', [1 0.5 0],'LineStyle', 'none');
    e2 = errorbar(min_sor_raftimedata,min_sor_rafdatanums,min_sor_rafdataerror,'LineStyle','none');
    e2.Color = 'black';
    e2.LineWidth = 1.4;
    p41 = plot(linspace(0,60),plus_sor_rafoutput(:,k),'Color',[0 0.5 1],'LineWidth', 2);
    p41.Color(4) = 0.25;
    p42 = plot(linspace(0,60),plus_sor_raf_optimalfit,'Color',[0 0.5 1],'LineWidth', 2);
    p43 = plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor', [0 0.5 1], 'LineStyle', 'none');
    e3 = errorbar(plus_sor_raftimedata,plus_sor_rafdatanums,plus_sor_rafdataerror,'LineStyle','none');
    e3.Color = 'black';
    e3.LineWidth = 1.4;
    legend([p32 p42],{'-sor membrane Raf1' '+sor membrane Raf1'});
    xlabel('Time (min)');
    ylabel('Concentration (molec/cell)');

    set(0, 'CurrentFigure', f3);
    p51 = plot(linspace(0,60),min_sor_perkoutput(:,k),'Color',[0.2 0.8 0.8],'LineWidth', 2);
    hold on
    p51.Color(4) = 0.25;
    p52 = plot(linspace(0,60),min_sor_perk_optimalfit,'Color',[0.2 0.8 0.8],'LineWidth', 2);
    p53 = plot(min_sor_perktimedata,min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor',[0.2 0.8 0.8],'MarkerEdgeColor', [0.2 0.8 0.8],'LineStyle', 'none');
    e4 = errorbar(min_sor_perktimedata,min_sor_perkdatanums,min_sor_perkdataerror,'LineStyle','none');
    e4.Color = 'black';
    e4.LineWidth = 1.4;
    
    p61 = plot(linspace(0,60),min_sor_pmekoutput(:,k),'Color',[236/255 0 140/255],'LineWidth', 2);
    p61.Color(4) = 0.25;
    p62 = plot(linspace(0,60),min_sor_pmek_optimalfit,'Color',[236/255 0 140/255],'LineWidth', 2);
    p63 = plot(min_sor_pmektimedata,min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor',[236/255 0 140/255],'MarkerEdgeColor', [236/255 0 140/255],'LineStyle', 'none');
    e5 = errorbar(min_sor_pmektimedata,min_sor_pmekdatanums,min_sor_pmekdataerror,'LineStyle','none');
    e5.Color = 'black';
    e5.LineWidth = 1.4;
    legend([p52 p62],{'-sor pERK' '-sor pMEK'});
    xlabel('Time (min)');
    ylabel('Concentration (molec/cell)');
    
end
hold off
%% Plot histogram of log-scaled fitted parameter values
figure
for m = 1:length(wanted_param2)
    subplot(8,5,m)
    histogram(log(fittedparams(:,m))); 
    xlabel(['log(' names{m} ')']);
end

%spider plot of fitted parameters
PLSR_cat_wanted = cell(length(wanted_param2),1)
for m = 1:length(wanted_param2)
    PLSR_cat_wanted{m} = names_all{wanted_param2(m)};
end
s = spider_plot(fittedparams,'FillOption', 'on','AxesLabels', PLSR_cat_wanted);

%% Parameter estimation (single fit) using a subset of parameters with VIP scores > 1  -> Nelder-Mead
objfunc_single = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata',wanted_param2) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata',wanted_param2) - min_sor_rafdatanums).^2./ (min_sor_rafdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param2) - plus_sor_rasdatanums).^2./ (plus_sor_rasdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param2) - plus_sor_rafdatanums).^2./ (plus_sor_rafdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata',wanted_param2) - min_sor_perkdatanums).^2./ (min_sor_perkdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param2) - min_sor_pmekdatanums).^2./ (min_sor_pmekdataerror).^2));

LB_single = zeros(1,length(wanted_param2))';
tic
[fittedparams_single,fval,exitflag,output]  = fminsearchbnd(objfunc_single,my_wanted_params,LB_single,[],options);
toc

[min_sor_rasoutput_single,T, min_sor_ras_allValues_single, param1_single, param_vary1_single] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'min_sor_rastimedata',wanted_param2);
[plus_sor_rasoutput_single,~, plus_sor_ras_allValues_single, param2_single, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'plus_sor_rastimedata',wanted_param2);
[min_sor_rafoutput_single,~, min_sor_raf_allValues_single, param3_single, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'min_sor_raftimedata',wanted_param2);
[plus_sor_rafoutput_single,~, plus_sor_raf_allValues_single, param4_single, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'plus_sor_raftimedata',wanted_param2);
[min_sor_perkoutput_single,~, min_sor_perk_allValues_single, param5_single, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'min_sor_perktimedata',wanted_param2);
[min_sor_pmekoutput_single,~, min_sor_pmek_allValues_single, param6_single, ~] = fullEGFR9_onemodel_fit_40_v2(fittedparams_single, linspace(0,60),'min_sor_pmektimedata',wanted_param2);


figure
plot(min_sor_rastimedata,min_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor',  [1 0.2 0.2],'LineStyle', 'none');
hold on
e1 = errorbar(min_sor_rastimedata,min_sor_rasdatanums,min_sor_rasdataerror,'LineStyle','none');
e1.Color = 'black';
e1.LineWidth = 1.4;
plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
p12 = plot(linspace(0,60),min_sor_rasoutput_single,'Color', [1 0.2 0.2],'LineWidth',2);
p22 = plot(linspace(0,60), plus_sor_rasoutput_single,'Color',[0 0.6 0.3],'LineWidth',2);
ax = gca;
ax.FontSize = 14; 
xlabel('time (min)','FontSize',14);
ylabel('Concentration (molec/cell)','FontSize',14);
legend([p12 p22],{'-sor GTP-bound Ras' '+sor GTP-bound Ras'},'FontSize',14);
hold off
saveas(gcf,'FigS2_subset_fit_RASGTP.pdf')

figure
plot(min_sor_raftimedata,min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor',[1 0.5 0],'LineStyle', 'none');
hold on
e2 = errorbar(min_sor_raftimedata,min_sor_rafdatanums,min_sor_rafdataerror,'LineStyle','none');
e2.Color = 'black';
e2.LineWidth = 1.4;
plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor', [0 0.5 1],'LineStyle', 'none');
e3 = errorbar(plus_sor_raftimedata,plus_sor_rafdatanums,plus_sor_rafdataerror,'LineStyle','none');
e3.Color = 'black';
e3.LineWidth = 1.4;
p32 = plot(linspace(0,60), min_sor_rafoutput_single,'Color',[1 0.5 0],'LineWidth',2);
p42 = plot(linspace(0,60), plus_sor_rafoutput_single,'Color',[0 0.5 1],'LineWidth',2);
ax = gca;
ax.FontSize = 14; 
xlabel('time (min)','FontSize',14);
ylabel('Concentration (molec/cell)','FontSize',14);
legend([p32 p42],{'-sor membrane Raf1' '+sor membrane Raf1'},'FontSize',14);
hold off
saveas(gcf,'FigS2_subset_fit_membRAF1.pdf')

figure
plot(min_sor_perktimedata,min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor', [0.2 0.8 0.8], 'MarkerEdgeColor', [0.2 0.8 0.8],'LineStyle', 'none');
hold on
e4 = errorbar(min_sor_perktimedata,min_sor_perkdatanums,min_sor_perkdataerror,'LineStyle','none');
e4.Color = 'black';
e4.LineWidth = 1.4;
plot(min_sor_pmektimedata,min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor', [236/255 0 140/255], 'MarkerEdgeColor', [236/255 0 140/255],'LineStyle', 'none');
e5 = errorbar(min_sor_pmektimedata,min_sor_pmekdatanums,min_sor_pmekdataerror,'LineStyle','none');
e5.Color = 'black';
e5.LineWidth = 1.4;
p52 = plot(linspace(0,60),min_sor_perkoutput_single,'Color',[0.2 0.8 0.8],'LineWidth',2);
p62 = plot(linspace(0,60), min_sor_pmekoutput_single,'Color',[236/255 0 140/255],'LineWidth',2);
ax = gca;
ax.FontSize = 14; 
xlabel('time (min)','FontSize',14)
ylabel('Concentration (molec/cell)','FontSize',14);
legend([p52 p62],{'-sor pERK' '-sor pMEK'},'FontSize',14);       
hold off
saveas(gcf,'FigS2_subset_fit_MEK_ERK.pdf')

%% Parameter estimation (single fit) using a subset of parameters with VIP scores > 1  -> Nelder-Mead, RAF1 params excluded

%% Parameter estimation (single fit) using a subset of parameters with VIP scores > 1 -> particleswarm
rng(1)
paramlist; % load best-fit values from SloppyCell
objfunc_single_all = @(input) obj_fnc_pswarm_subset(wanted_param2, input);

my_wanted_params = params(wanted_param2)';
LB_single_all = my_wanted_params .* 0.001;

%LB_single_all = zeros(1,length(wanted_param2))';
%LB_single_all = zeros(1,length(wanted_param2));

%LB_single_all = min(fittedparams); %get min value for each param

%UB_single_all = my_wanted_params .* 100; %zeros(1,length(wanted_param))' + 10^6;
UB_single_all = [];
%UB_single_all = max(fittedparams); %get min value for each param


swarmsize = 50;
%options = optimoptions('particleswarm','SwarmSize', swarmsize, 'UseParallel', false, 'PlotFcn','pswplotbestf', 'FunctionTolerance', 1e-4);
%options = optimoptions('particleswarm','SwarmSize', 400, 'HybridFcn',@fmincon, 'PlotFcn','pswplotbestf');
options = optimoptions('particleswarm','SwarmSize', swarmsize, 'UseParallel', false, 'PlotFcn',@pswplotbestf, 'OutputFcn',@pswoutfun, 'FunctionTolerance', 1e-2,'MaxIterations',1000);


%nvars = length(wanted_param2);
nvars = length(wanted_param2);


tic
%fittedparams = particleswarm(objfunc_single_all, nvars, LB_single_all, UB_single_all, options);
[fittedparams_pswarm, fval, exitflag, pswarm_output] = particleswarm(objfunc_single_all, nvars, LB_single_all, UB_single_all, options);
plot_est(fittedparams_pswarm, wanted_param2, 'pswarm', swarmsize)
cost = zeros(length(swarm_vals),1);
for i = 1:length(swarm_vals)
    
    for k=1:2

        cost(i,k) = 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), min_sor_rastimedata,'min_sor_rastimedata',wanted_param) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
                         + sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), min_sor_raftimedata,'min_sor_raftimedata',wanted_param) - min_sor_rafdatanums).^2 ./ (min_sor_rafdataerror).^2)...
                         + sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param) - plus_sor_rasdatanums).^2 ./ (plus_sor_rasdataerror).^2) ...
                         + sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param) - plus_sor_rafdatanums).^2 ./ (plus_sor_rafdataerror).^2) ...
                         + sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), min_sor_perktimedata,'min_sor_perktimedata',wanted_param) - min_sor_perkdatanums).^2 ./ (min_sor_perkdataerror).^2) ...
                         + sum((fullEGFR9_onemodel_fit_40_v2(swarm_vals{i,2}(k,:), min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param) - min_sor_pmekdatanums).^2 ./ (min_sor_pmekdataerror).^2));
    end
        
end


toc


%% Single best fit using all 40 parameters
my_wanted_params = params(wanted_param)';
objfunc_single_all = @(x) 1/2 * (sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_rastimedata,'min_sor_rastimedata',wanted_param) - min_sor_rasdatanums).^2 ./(min_sor_rasdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_raftimedata,'min_sor_raftimedata',wanted_param) - min_sor_rafdatanums).^2./ (min_sor_rafdataerror).^2)...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_rastimedata,'plus_sor_rastimedata',wanted_param) - plus_sor_rasdatanums).^2./ (plus_sor_rasdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, plus_sor_raftimedata,'plus_sor_raftimedata',wanted_param) - plus_sor_rafdatanums).^2./ (plus_sor_rafdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_perktimedata,'min_sor_perktimedata',wanted_param) - min_sor_perkdatanums).^2./ (min_sor_perkdataerror).^2) ...
    + sum((fullEGFR9_onemodel_fit_40_v2(x, min_sor_pmektimedata,'min_sor_pmektimedata',wanted_param) - min_sor_pmekdatanums).^2./ (min_sor_pmekdataerror).^2));

LB_single_all = zeros(1,length(wanted_param))';
tic
[fittedparams_single_all,fval,exitflag,output]  = fminsearchbnd(objfunc_single_all,my_wanted_params,LB_single_all,[],options);
toc


%% Plot different rates of convergence of objective function values according to different numbers of parameters selected by PLSR SA 


save estimation_main.mat