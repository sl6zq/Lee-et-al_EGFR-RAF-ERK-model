%% Comparison of univariate sensitivity analysis vs. PLSR SA (PLSR coefficients) for individual model outputs (RAS-GTP, memb. RAF1, pMEK, pERK) 
clc
clear all
% setting to save ODE solutions for dependent variable in PLSR SA 
savedata = false;
%% Univariate sensitivity analysis
% Load best-fit paramter values
paramlist;

paramnames = {'k_{E,f}'; 'k_{catE}'; 'k_{pMEK}'; 'k_{nfpSOS}'; 'k_{B,r}'; 'k_{B,f}'; 'k_{dpERK}'; 'k_{R1,r}'; 'k_{R1,f}'; 'k_{nfpBR1}'; 'k_{pR1}'; 'k_{SOS,r}'; 'k_{SOS,f}'; 'k_{Rgneslow}'; 'k_{Son}'; 'k_{G2SOS,r}'; 'k_{G2SOS,f}'; 'k_{iR1,r}'; 'k_{iR1,f}'; 'k_{Soff}'; 'k_{pERK}'; 'k_{Rhydro}';
'K_{mgneslow}'; 'k_{dpMEK}'; 'k_{nfpiR1r}'; 'k_{dR1,r}'; 'k_{iB,r}'; 'k_{dR1,f}'; 'k_{iB,f}'; 'k_{fpB,r}'; 'k_{dpSOS}';
'k_{dp}'; 'k_{dE,r}'; 'k_{dE,f}'; 'k_{dpR1}'; 'k_{G2,r}'; 'k_{G2,f}'; 'k_{nfpiB,r}'; 'k_{fpR1,r}'; 'k_{E,r}'};

speciesindex = [184;157;21;45];
speciesnames = {'GTP-RAS'; 'Membrane-bound RAF1'; 'pMEK'; 'pERK'};

labels = paramnames';

timeSpan = linspace(0.0, 60.0);

[T,~,~,~,~,allValues_minsor] = fullEGFR9_onemodel(timeSpan, yinit, params, 'min_sor', '');
[~,~,~,~,~,allValues_plussor] = fullEGFR9_onemodel(timeSpan, yinit, params, 'plus_sor', '');

baseYave_minsor = zeros(length(speciesindex),1);
baseYave_plussor = zeros(length(speciesindex),1);
baseYave_2 = zeros(length(speciesindex),1);
tryout = zeros(length(T),length(speciesindex));
for i = 1:length(speciesindex)
    tryout(:,i) = allValues_minsor(:,speciesindex(i));
    baseYave_2(i) = trapz(allValues_minsor(:,speciesindex(i)),T);
    baseYave_2(i) = trapz(allValues_plussor(:,speciesindex(i)),T);
end


baseYave_minsor = baseYave_minsor./(T(end)-T(1));
baseYave_plussor = baseYave_plussor./(T(end)-T(1));
baseYave_2 = baseYave_2./(T(end)-T(1));
baseYmax_minsor = max(allValues_minsor(:,speciesindex))';
baseYmax_plussor = max(allValues_plussor(:,speciesindex))';

cell_wanted = (num2cell(wanted_param));

names_nums = horzcat(cell_wanted, paramnames);
store = zeros(length(T),length(wanted_param),2);
inc_param_minsor = params;
dec_param_minsor = params;
inc_param_plussor = params;
dec_param_plussor = params;
aveYsens = zeros(length(speciesindex), length(wanted_param));
maxYsens = zeros(length(speciesindex), length(wanted_param));
for i = 1: length(wanted_param)
    dec_param_minsor(wanted_param(i)) = 0.05 .* params(wanted_param(i));
    inc_param_minsor(wanted_param(i)) = 5 .* params(wanted_param(i));
    
    dec_param_plussor(wanted_param(i)) = 0.05 .* params(wanted_param(i));
    inc_param_plussor(wanted_param(i)) = 5 .* params(wanted_param(i));
   
    %DECREASE
    [decT_minsor,~,~,~,~,allValuesdec_minsor] = fullEGFR9_onemodel(T, yinit, dec_param_minsor, 'min_sor','');
    %INCREASE
    [incT_minsor,~,~,~,~,allValuesinc_minsor] = fullEGFR9_onemodel(T, yinit, inc_param_minsor, 'min_sor','');
    %DECREASE
    [decT_plussor,~,~,~,~,allValuesdec_plussor] = fullEGFR9_onemodel(T, yinit, dec_param_plussor, 'plus_sor','');
    %INCREASE
    [incT_plussor,~,~,~,~,allValuesinc_plussor] = fullEGFR9_onemodel(T, yinit, inc_param_plussor, 'plus_sor','');
    
    %Y SPECIES CALCS
    %calculate average value of each species over simulation time
    for j = 1:length(speciesindex)
        decYave_minsor(j,1) = trapz(allValuesdec_minsor(:,speciesindex(j)),decT_minsor);
        incYave_minsor(j,1) = trapz(allValuesinc_minsor(:,speciesindex(j)),incT_minsor);
        decYave_plussor(j,1) = trapz(allValuesdec_plussor(:,speciesindex(j)),decT_plussor);
        incYave_plussor(j,1) = trapz(allValuesinc_plussor(:,speciesindex(j)),incT_plussor);
       
    end
    decYave_minsor = decYave_minsor./(decT_minsor(end)-decT_minsor(1));
    incYave_minsor = incYave_minsor./(incT_minsor(end)-incT_minsor(1));
    decYave_plussor = decYave_plussor./(decT_plussor(end)-decT_plussor(1));
    incYave_plussor = incYave_plussor./(incT_plussor(end)-incT_plussor(1));
    
    decYmax_minsor = max(allValuesdec_minsor(:,speciesindex))';
    incYmax_minsor = max(allValuesinc_minsor(:,speciesindex))';
    decYmax_plussor = max(allValuesdec_plussor(:,speciesindex))';
    incYmax_plussor = max(allValuesinc_plussor(:,speciesindex))';
    
    %calculate order-of-magnitude differences of species averages compared
    %to base model case
    decYavediff_minsor = log10(decYave_minsor./baseYave_minsor);
    incYavediff_minsor = log10(incYave_minsor./baseYave_minsor);
    decYavediff_plussor = log10(decYave_plussor./baseYave_plussor);
    incYavediff_plussor = log10(incYave_plussor./baseYave_plussor);

    decYmaxdiff_minsor = log10(decYmax_minsor./baseYmax_minsor);
    incYmaxdiff_minsor = log10(incYmax_minsor./baseYmax_minsor);
    decYmaxdiff_plussor = log10(decYmax_plussor./baseYmax_plussor);
    incYmaxdiff_plussor = log10(incYmax_plussor./baseYmax_plussor);
    
    aveYsens_minsor(:,i) = abs(decYavediff_minsor)+abs(incYavediff_minsor);
    maxYsens_minsor(:,i) = abs(decYmaxdiff_minsor)+abs(incYmaxdiff_minsor);
    aveYsens_plussor(:,i) = abs(decYavediff_plussor)+abs(incYavediff_plussor);
    maxYsens_plussor(:,i) = abs(decYmaxdiff_plussor)+abs(incYmaxdiff_plussor);
    
    %Reset model parameters
    inc_param_minsor = params;
    dec_param_minsor = params;
    inc_param_plussor = params;
    dec_param_plussor = params;

end

maxaveYsens_minsor = max(aveYsens_minsor,[],2);
maxmaxYsens_minsor = max(maxYsens_minsor,[],2);
maxaveYsens_plussor = max(aveYsens_plussor,[],2);
maxmaxYsens_plussor = max(maxYsens_plussor,[],2);
for i = 1:length(maxaveYsens_minsor)
    aveYsens_minsor(i,:) = aveYsens_minsor(i,:)./maxaveYsens_minsor(i).*100;
    maxYsens_minsor(i,:) = maxYsens_minsor(i,:)./maxmaxYsens_minsor(i).*100;
    aveYsens_plussor(i,:) = aveYsens_plussor(i,:)./maxaveYsens_plussor(i).*100;
    maxYsens_plussor(i,:) = maxYsens_plussor(i,:)./maxmaxYsens_plussor(i).*100;
end 


figure
subplot(2,1,1);
hm2 = heatmap(labels,speciesnames, maxYsens_minsor, 'FontSize', 11.5);
colormap(jet);
ax = gca;
title('-sorafenib');
hm2.FontSize = 12;
subplot(2,1,2)
hm3 = heatmap(labels, speciesnames, maxYsens_plussor, 'FontSize', 11.5);
colormap(jet);
ax = gca;
title('+sorafenib');
hm3.FontSize = 12;
set(gcf, 'color','w');
set(gcf,'Position',[10, 10, 1600, 410]);

%% Define parameters
rng(0);                      %random seed for SA
paramlist;                   %baseline parameters (best fit param values from SloppyCell)
sampling         = 'random'; %sampling method for SA, choices: random, LHS, log-random, log-LHS
% ODE settings
num_Param        = 3000;     % at least 2 for cv partition, num_Param >= 122
timeSpan         = 0:1:60;   % simulate first 60 min after EGF treatment
ncomp            = 40;       % PLSR model with 40 componenets (max. possible # of components)

%% Get ODE solutions with 3000 parameter sets 
paramIn         = cell(1,6);
paramOut        = cell(1,12);

[output] = evaluateparam_one(wanted_param, params, timeSpan, yinit, num_Param, sampling);

if savedata
    save ODE_slns_3000.mat
end
%% 1) PLSR model
% PLSR w/ 6 maximum outputs: +/-sor RASGTP, +/-memb RAF1, pMEK, pERK 
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(num_Param, sampling, output.my_other_params_plussor, output.my_y_max_total, ncomp, wanted_param);

%% Compare PLSR coefficients with univariate SA (Fig. S3)
timeSpan                  = 0:1:60;
coeff_percentmax          = zeros(40,1);
for i=1:6
    coeff_percentmax(:,i) = abs(BETA(2:end,i)'/max(abs(BETA(2:end,i)))) * 100;
end
yvalues                   = PLSR_cat;
min_sor_speciesnames      = {'GTP-RAS'; 'Membrane-bound RAF1';'pMEK'; 'pERK'};

sens_all                  = cat(2,maxYsens_minsor', coeff_percentmax(:,3:6));
all_speciesnames          = {'Uni GTP-loaded RAS'; 'Uni Membrane-bound RAF1';'Uni pMEK'; 'Uni pERK'; 'Multi GTP-loaded RAS'; 'Multi Membrane-bound RAF1';'Multi pMEK'; 'Multi pERK'};

figure
h1 = heatmap(yvalues', all_speciesnames, sens_all','CellLabelColor','none','Colormap',jet,'FontSize',9);
ax = gca;
set(gca,'FontSize',9);