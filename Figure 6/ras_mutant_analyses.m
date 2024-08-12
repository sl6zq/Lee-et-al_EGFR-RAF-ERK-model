clear all
clc
% If savedata is true, save all plots and data in .mat files
savedata = false;
% 1) 6A: Signaling dynamics at key pathway nodes for WT RAS vs. mutant RAS
% 2) 6B: PLSR SA VIP scores with difference in time-integrated pERK levels
% between WT vs. mutant
% 3) 6C: pERK dynamics for WT RAS, mutant RAS, WT RAS with RAF1 * 10, mutant RAS with RAF1 * 10
% 3) 6D: pERK dynamics no feedback
% 4) 6E: Difference in time-integrated levels of pERK (WT vs. mutant for baseline (low) [RAF1] vs. [RAF1] * 10)
% 5) 6F: Difference in time-integrated levels of all model species (WT vs. mutant for baseline (low) [RAF1] vs. [RAF1] * 10)
timeSpan = 0:1:60;
%% 6A
traj_sormin    = load('ens-trajectories-sorminus');
traj_sorplus   = load('ens-trajectories-sorplus');
mutant_sormin  = load('ens-mutant-trajectories-sorminus');
mutant_sorplus = load('ens-mutant-trajectories-sorplus');

diff_RasGTP_mu_minsor  = round(max(mutant_sormin.RasGTP_mu) / max(traj_sormin.RasGTP_mu),2);
diff_Raf1_mu_minsor    = round(max(mutant_sormin.Raf1_mu) / max(traj_sormin.Raf1_mu),2);
diff_pMEK_mu_minsor    = round(max(mutant_sormin.pMEK_mu) / max(traj_sormin.pMEK_mu),2);
diff_pERK_mu_minsor    = round(max(mutant_sormin.pERK_mu) / max(traj_sormin.pERK_mu),2);

diff_RasGTP_mu_plussor = round(max(mutant_sorplus.RasGTP_mu) / max(traj_sorplus.RasGTP_mu),2);
diff_Raf1_mu_plussor   = round(max(mutant_sorplus.Raf1_mu) / max(traj_sorplus.Raf1_mu),2);

figure
m1 = plot(traj_sormin.t,mutant_sormin.RasGTP_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b1 = plot(traj_sormin.t,mutant_sormin.RasGTP_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p1 = patch([traj_sormin.t fliplr(traj_sormin.t)], [mutant_sormin.RasGTP_mu-mutant_sormin.RasGTP_sigma  fliplr(mutant_sormin.RasGTP_mu + mutant_sormin.RasGTP_sigma)], [248/255 118/255 109/255], 'FaceAlpha',0.3);
p1.EdgeColor = 'none'; title('RAS-GTP', 'FontSize', 15);
m2 = plot(traj_sormin.t,traj_sormin.RasGTP_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
p2 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.RasGTP_mu-traj_sormin.RasGTP_sigma  fliplr(traj_sormin.RasGTP_mu + traj_sormin.RasGTP_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3);
b2= plot(traj_sormin.t,traj_sormin.RasGTP_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p2.EdgeColor = 'none'; text(max(traj_sormin.t)*0.7, max([mutant_sormin.RasGTP_best traj_sormin.RasGTP_best]), ['\Delta Max.: ' num2str(diff_RasGTP_mu_minsor)]);
xlabel('Time (min)'); ylabel('Concentration (molec/cell)');
legend([b1 b2],{'KRAS G12V -sor GTP-bound Ras ' 'WT -sor GTP-bound Ras '},'FontSize',9,'Location', 'northoutside');

ax = gca;
ax.YAxis.Exponent = 0;
hold off

figure
m3 = plot(traj_sormin.t,mutant_sormin.Raf1_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b3 = plot(traj_sormin.t,mutant_sormin.Raf1_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p3 = patch([traj_sormin.t fliplr(traj_sormin.t)], [mutant_sormin.Raf1_mu-mutant_sormin.Raf1_sigma  fliplr(mutant_sormin.Raf1_mu + mutant_sormin.Raf1_sigma)], [248/255 118/255 109/255], 'FaceAlpha',0.3);
p3.EdgeColor = 'none'; title('Membrane RAF1', 'FontSize', 15);
m4 = plot(traj_sormin.t,traj_sormin.Raf1_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
b4 = plot(traj_sormin.t,traj_sormin.Raf1_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p4 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.Raf1_mu-traj_sormin.Raf1_sigma  fliplr(traj_sormin.Raf1_mu + traj_sormin.Raf1_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3);
p4.EdgeColor = 'none'; xlabel('Time (min)'); ylabel('Concentration (molec/cell)'); 
text(max(traj_sormin.t)*0.7, max([mutant_sormin.Raf1_best traj_sormin.Raf1_best]), ['\Delta Max.: ' num2str(diff_Raf1_mu_minsor)]);
legend([b3 b4],{'KRAS G12V -sor membrane Raf1' 'WT -sor membrane Raf1'},'FontSize',9,'Location', 'northoutside');
hold off

figure
m5 = plot(traj_sormin.t,mutant_sormin.pMEK_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b5 = plot(traj_sormin.t,mutant_sormin.pMEK_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p5 = patch([traj_sormin.t fliplr(traj_sormin.t)], [mutant_sormin.pMEK_mu-mutant_sormin.pMEK_sigma  fliplr(mutant_sormin.pMEK_mu + mutant_sormin.pMEK_sigma)],[248/255 118/255 109/255], 'FaceAlpha',0.3)
p5.EdgeColor = 'none';
m6 = plot(traj_sormin.t,traj_sormin.pMEK_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
b6 = plot(traj_sormin.t,traj_sormin.pMEK_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p6 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.pMEK_mu-traj_sormin.pMEK_sigma  fliplr(traj_sormin.pMEK_mu + traj_sormin.pMEK_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3)
p6.EdgeColor = 'none';
xlabel('Time (min)'); ylabel('Concentration (molec/cell)');
ylim([-0.5*10^5 3.2*10^5]);title('pMEK', 'FontSize', 15);
text(max(traj_sormin.t)*0.7, max([mutant_sormin.pMEK_best traj_sormin.pMEK_best])*1.2, ['\Delta Max.: ' num2str(diff_pMEK_mu_minsor)]);
legend([b5 b6],{'KRAS G12V -sor pMEK' 'WT -sor pMEK'},'FontSize',9,'Location', 'northoutside');       
hold off

figure
m5 = plot(traj_sormin.t,mutant_sormin.pERK_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b5 = plot(traj_sormin.t,mutant_sormin.pERK_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p5 = patch([traj_sormin.t fliplr(traj_sormin.t)], [mutant_sormin.pERK_mu-mutant_sormin.pERK_sigma  fliplr(mutant_sormin.pERK_mu + mutant_sormin.pERK_sigma)],[248/255 118/255 109/255], 'FaceAlpha',0.3)
p5.EdgeColor = 'none';
m6 = plot(traj_sormin.t,traj_sormin.pERK_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
b6 = plot(traj_sormin.t,traj_sormin.pERK_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p6 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.pERK_mu-traj_sormin.pERK_sigma  fliplr(traj_sormin.pERK_mu + traj_sormin.pERK_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3)
p6.EdgeColor = 'none';
xlabel('Time (min)');
ylabel('Concentration (molec/cell)');
text(max(traj_sormin.t)*0.7, max([mutant_sormin.pERK_best traj_sormin.pERK_best])*1.1, ['\Delta Max.: ' num2str(diff_pERK_mu_minsor)]);
legend([b5 b6],{'KRAS G12V -sor pERK' 'WT -sor pERK'},'FontSize',9,'Location', 'northoutside');
ylim([-0.5*10^5 3.2*10^5]); title('pERK', 'FontSize', 15);
hold off

figure
m5 = plot(traj_sorplus.t,mutant_sorplus.Raf1_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b5 = plot(traj_sorplus.t,mutant_sorplus.Raf1_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p5 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [mutant_sorplus.Raf1_mu-mutant_sorplus.Raf1_sigma  fliplr(mutant_sorplus.Raf1_mu + mutant_sorplus.Raf1_sigma)],[248/255 118/255 109/255], 'FaceAlpha',0.3)
p5.EdgeColor = 'none';
m6 = plot(traj_sorplus.t,traj_sorplus.Raf1_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
b6 = plot(traj_sorplus.t,traj_sorplus.Raf1_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p6 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [traj_sorplus.Raf1_mu-traj_sorplus.Raf1_sigma  fliplr(traj_sorplus.Raf1_mu + traj_sorplus.Raf1_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3)
p6.EdgeColor = 'none';
xlabel('Time (min)'); ylabel('Concentration (molec/cell)');
ylim([0 12000]);ax = gca;ax.YAxis.Exponent = 0;
title('Membrane RAF1', 'FontSize', 15);
text(max(traj_sorplus.t)*0.7, max([mutant_sorplus.Raf1_best traj_sorplus.Raf1_best])*1.1, ['\Delta Max.: ' num2str(diff_Raf1_mu_plussor)]);
legend([b5 b6],{'KRAS G12V +sor membrane RAF1' 'WT +sor membrane RAF1'},'FontSize',9,'Location', 'northoutside');       
hold off

figure
m5 = plot(traj_sorplus.t,mutant_sorplus.RasGTP_mu,'--','color', [248/255 118/255 109/255], 'LineWidth', 2);
hold on
b5 = plot(traj_sorplus.t,mutant_sorplus.RasGTP_best,'-','color', [248/255 118/255 109/255], 'LineWidth', 2);
p5 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [mutant_sorplus.RasGTP_mu-mutant_sorplus.RasGTP_sigma  fliplr(mutant_sorplus.RasGTP_mu + mutant_sorplus.RasGTP_sigma)],[248/255 118/255 109/255], 'FaceAlpha',0.3)
p5.EdgeColor = 'none';
m6 = plot(traj_sorplus.t,traj_sorplus.RasGTP_mu,'--','color', [0/255 191/255 196/255], 'LineWidth', 2);
b6 = plot(traj_sorplus.t,traj_sorplus.RasGTP_best,'-','color', [0/255 191/255 196/255], 'LineWidth', 2);
p6 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [traj_sorplus.RasGTP_mu-traj_sorplus.RasGTP_sigma  fliplr(traj_sorplus.RasGTP_mu + traj_sorplus.RasGTP_sigma)], [0/255 191/255 196/255], 'FaceAlpha',0.3)
p6.EdgeColor = 'none';
xlabel('Time (min)');
ylabel('Concentration (molec/cell)');title('RAS-GTP', 'FontSize', 15);
text(max(traj_sorplus.t)*0.7, max([mutant_sorplus.RasGTP_best traj_sorplus.RasGTP_best])*1.1, ['\Delta Max.: ' num2str(diff_RasGTP_mu_plussor)]);
legend([b5 b6],{'KRAS G12V +sor GTP-bound RAS' 'WT +sor GTP-bound RAS'},'FontSize',9,'Location', 'northoutside');       
hold off

%% 6B: Parameter sensitivity analysis of fold change in pERK for KRAS mutant vs. WT RAS 
% Load baseline parameter values
paramlist;
paramlist_krasmutant; %16-fold decrease in kRhydro

%specify number of parameter sets 
number = 3000;
%specify parameter sampling method
sampling = 'random';
%set random seed
rng(1);
% Solve ODEs with 3000 parameter arrays for wild-type RAS
WT_RAS   = evaluateparam_one(wanted_param_wt, params, timeSpan, yinit, number, sampling);
% Solve ODEs with 3000 parameter arrays for KRAS mutant
KRASG12V = evaluateparam_one(wanted_param_mutant, params_mutant, timeSpan, yinit_mutant, number, sampling);
%Set change in pERK as Y in PLSR
pERK_change = abs(WT_RAS.my_y_max_total{6} - KRASG12V.my_y_max_total{6}); % change caused by parameter variations

ncomp = 10;
[plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(number, sampling, WT_RAS.my_other_params_minsor, pERK_change, ncomp, wanted_param_wt);

plot_plsr(PLSR_cat_mutant, number, sampling, Xloadings, vipScores, PCTVAR, 'vip', 'pERK difference, WT/KRASG12V');
%% 6C: Test the effect of increasing RAF1 expression in wild-type RAS & KRAS mutant 
timeSpan = 0:1:60;
%WT Ras w/o sorafenib
[T,~,~,params_minsor,allNames,allValues_wt]             = fullEGFR9_onemodel(timeSpan, yinit, params, 'min_sor', 'no');
%WT Ras w/ sorafenib
[~,~,~,params_plussor,~,~]                              = fullEGFR9_onemodel(timeSpan, yinit, params, 'plus_sor', 'no');
%mutant w/o sorafenib
[~,~,~,params_mutant,~,allValues_mutant]                = fullEGFR9_onemodel(timeSpan, yinit_mutant, params_mutant,'min_sor', 'no');

yinit_raf_highexp     = yinit_mutant;
yinit_raf_highexp(15) = yinit_mutant(15)*10; % RAF1 increase 10-fold

[~,~,~,params_raf_highexp,~,allValues_raf_highexp] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_mutant,'min_sor', 'no'); %KRAS G12V with 10-fold increase in [RAF1]
[~,~,~,params_wt_raf_highexp,~,allValues_wt_raf_highexp] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params,'min_sor', 'no'); %wild-type RAS with 10-fold increase in [RAF1]

h1 = figure;
set(h1,'Position',[100 100 400 300]);
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
if savedata == 1
    savefig('mutant_wt_raf_init_change.fig');
end
fold_diff1    = abs(trapz(T,allValues_wt(:,45)) - trapz(T,allValues_mutant(:,45)));
fold_diff2    = abs(trapz(T,allValues_wt_raf_highexp(:,45)) - trapz(T,allValues_raf_highexp(:,45)));
rafinc_effect = fold_diff2/fold_diff1;
%% 6D: Test feedback role (nfpSOS, nfpBR1)  
params_mutant_nofdback      = params_mutant;
params_mutant_nofdback(10)  = 0; %'knfpSOS'
params_mutant_nofdback(45)  = 0; %'knfpBR1'
params_mutant_nofdback(111) = 0; %'kfpBr'
params_mutant_nofdback(63)  = 0; %'kfpR1r'

params_nofdback             = params;
params_nofdback(10)         = 0; %'knfpSOS'
params_nofdback(45)         = 0; %'knfpBR1'
params_nofdback(111)        = 0; %'kfpBr'
params_nofdback(63)         = 0; %'kfpR1r'

[~,~,~,params_wt_raf_nofdbk,~,allValues_wt_RAS_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit, params_nofdback,'min_sor', 'no'); %WT RAS, no fdback
[~,~,~,params_raf_nofdbk,~,allValues_KRASG12V_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit, params_mutant_nofdback,'min_sor', 'no'); %KRASG12V, no fdback
[~,~,~,params_raf_highexp_nofdbk,~,allValues_raf_highexp_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_mutant_nofdback,'min_sor', 'no'); %KRAS G12V with 10-fold increase in [RAF1]
[~,~,~,params_wt_raf_highexp_nofdbk,~,allValues_wt_raf_highexp_nofdbk] = fullEGFR9_onemodel(timeSpan, yinit_raf_highexp, params_nofdback,'min_sor', 'no'); %wild-type RAS with 10-fold increase in [RAF1]


h1 = figure;
set(h1,'Position',[100 100 400 300]);
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
if savedata == 1
    savefig('mutant_wt_raf_init_change.fig');
end
fold_diff1_nofdbck    = abs(trapz(T,allValues_wt(:,45)) - trapz(T,allValues_mutant(:,45)));
fold_diff2_nofdbck    = abs(trapz(T,allValues_wt_raf_highexp_nofdbk(:,45)) - trapz(T,allValues_raf_highexp_nofdbk(:,45)));
rafinc_effect_nofdbck = fold_diff2_nofdbck/fold_diff1_nofdbck;
fold_diff2_nofdbck2   = abs(trapz(T,allValues_wt_RAS_nofdbk(:,45)) - trapz(T,allValues_KRASG12V_nofdbk(:,45)));
%% 6E
bar_x = categorical({'12,000', '120,000', '12,000 fdbck off', '120,000 fdbck off'});
bar_x = reordercats(bar_x, {'12,000', '120,000', '12,000 fdbck off', '120,000 fdbck off'});
bar_y = [fold_diff1 fold_diff2 fold_diff2_nofdbck2 fold_diff2_nofdbck];
figure
set(gcf,'Position',[100 100 300 400])
bar(bar_x, bar_y, 'BarWidth',0.6);
ylim([0 max(bar_y)*1.2]);
if savedata == 1
    savefig('mutant_wt_raf_tint_comparisons.fig');
end
%% Check concentrations for MAPK species
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
if savedata == 1
    figure
    tiledlayout(4,6, 'Padding', 'none', 'TileSpacing', 'compact');
    for i=1:length(pway_order_num)/2
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
end
%min-max scaling 
allValues_mutant_norm                                   = zeros(size(allValues_mutant,1), size(allValues_mutant,2)); 
allValues_wt_norm                                       = zeros(size(allValues_wt,1), size(allValues_wt,2)); 
for i=1:length(species_list)
    allValues_mutant_norm(:,i)                          = (allValues_mutant(:,species_list(i)) - min(allValues_mutant(:,species_list(i)))) ./ (max(allValues_mutant(:,species_list(i))) - min(allValues_mutant(:,species_list(i))));
    allValues_wt_norm(:,i)                              = (allValues_wt(:,species_list(i)) - min(allValues_wt(:,species_list(i)))) ./ (max(allValues_wt(:,species_list(i))) - min(allValues_wt(:,species_list(i))));
end

delta_values = zeros(length(species_list),1);
for i=1:length(species_list)
    delta_values(i)                                     = trapz(allValues_mutant(:,species_list(i))) ./ trapz(allValues_wt(:,species_list(i)));
end
delta_group                                             = [delta_values string(species_names)];
sorted_delta_group                                      = sortrows(delta_group,1,'descend');
sorted_delta_numbers                                    = str2double(sorted_delta_group(:,1));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
pway_delta_values                                       = zeros(length(pway_order_num),1);
for i=1:length(pway_order)
    pway_delta_values(i)                                = abs(trapz(allValues_mutant_norm(:,pway_order_num(i))) - trapz(allValues_wt_norm(:,pway_order_num(i))));
end

high_delta_index                                        = find(pway_delta_values > 10);   
[delta_group_pway_nonzero_idx, ~]                       = find(pway_delta_values > 0);  

figure
bar(pway_delta_values(delta_group_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
plot(xlim,[10 10],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold'); ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_pway_nonzero_idx),'xTickLabel',pway_order_group(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_pway_nonzero_idx),diag(pway_delta_values(delta_group_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [230/255 159/255 0/255]; %227, 111, 71
end
set(gcf,'Position',[100 100 750 500]);
if savedata ==1
    saveas(gcf,[pwd '/Plots/timeint_changes.pdf']);
end

figure
bar(pway_delta_values(delta_group_pway_nonzero_idx));
hold on
xlim=get(gca,'xlim');
plot(xlim,[10 10],'-.k','LineWidth',1);
xlabel('Node','FontWeight','bold');
ylabel({'Abs. diff. between time-integrated conc.'; 'in mutatnt vs. WT'},'FontWeight','bold');
set(gca, 'xTick', 1:length(delta_group_pway_nonzero_idx),'xTickLabel',pway_order_group(:,2),'XTickLabelRotation',45);
bh = bar(1:numel(delta_group_pway_nonzero_idx),diag(pway_delta_values(delta_group_pway_nonzero_idx)),'stacked','FaceColor', [153/255 153/255 255/255]);  %0, 155, 250
for k=1:length(high_delta_index)
    bh(high_delta_index(k)).FaceColor = [153/255 153/255 255/255]; %227, 111, 71
end
set(gcf,'Position',[100 100 750 500]);

%% with 10-fold increase in low abundance species' expression:
lowabund_species = 'RAF1';
if strcmp(lowabund_species, 'RAF1')
    yinit_raf_highexp                                                   = yinit_mutant;
    yinit_raf_highexp(15)                                               = yinit_mutant(15)*10; % RAF1 increase 10-fold
    params_raf_highexp                                                  = params;
    params_raf_highexp(115)                                             = params_raf_highexp(115) * 10; % RAF1 increase 10-fold
    params_mutant_raf_highexp                                           = params_mutant;
    params_mutant_raf_highexp(115)                                      = params_mutant_raf_highexp(115) * 10; % RAF1 increase 10-fold
elseif strcmp(lowabund_species, 'BRAF')
    yinit_raf_highexp                                                   = yinit_mutant;
    yinit_raf_highexp(33)                                               = yinit_mutant(33)*10; % BRAF increase 10-fold
    params_raf_highexp                                                  = params;
    params_raf_highexp(131)                                             = params_raf_highexp(131) * 10;
    params_mutant_raf_highexp                                           = params_mutant;
    params_mutant_raf_highexp(131)                                      = params_mutant_raf_highexp(131) * 10; % BRAF increase 10-fold
elseif strcmp(lowabund_species,'SOS')
    yinit_raf_highexp                                                   = yinit_mutant;
    yinit_raf_highexp(34)                                               = yinit_mutant(34)*10; % SOS increase 10-fold
    params_raf_highexp                                                  = params;
    params_raf_highexp(11)                                              = params_raf_highexp(11) * 10; % SOS increase 10-fold
    params_mutant_raf_highexp                                           = params_mutant;
    params_mutant_raf_highexp(11)                                       = params_mutant_raf_highexp(11) * 10; % SOS increase 10-fold
end

% WT-Ras dynamics
[~,~,~,out_params_wt_raf_highexp,~,allValues_wt_raf_highexp]        = fullEGFR9_onemodel(timeSpan,yinit_raf_highexp, params_raf_highexp,'min_sor','no'); %wild-type RAS with 10-fold increase in [RAF1]

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

%min-max scaling 
allValues_rafinc_mutant_norm                                        = zeros(size(allValues_mut_raf_highexp,1), size(allValues_mut_raf_highexp,2)); 
allValues_rafinc_wt_norm                                            = zeros(size(allValues_wt_raf_highexp,1), size(allValues_wt_raf_highexp,2)); 
for i=1:length(species_list)
    allValues_rafinc_mutant_norm(:,i)                               = (allValues_mut_raf_highexp(:,species_list(i)) - min(allValues_mut_raf_highexp(:,species_list(i)))) ./ (max(allValues_mut_raf_highexp(:,species_list(i))) - min(allValues_mut_raf_highexp(:,species_list(i))));
    allValues_rafinc_wt_norm(:,i)                                   = (allValues_wt_raf_highexp(:,species_list(i)) - min(allValues_wt_raf_highexp(:,species_list(i)))) ./ (max(allValues_wt_raf_highexp(:,species_list(i))) - min(allValues_wt_raf_highexp(:,species_list(i))));
end

% Plot in descending order
delta_values_rafinc                                                 = zeros(length(species_list),1);
for i=1:length(species_list)
    delta_values_rafinc(i)                                          = abs(trapz(allValues_rafinc_mutant_norm(:,species_list(i))) - trapz(allValues_rafinc_wt_norm(:,species_list(i))));
end
delta_group_rafinc                                                  = [delta_values_rafinc string(species_names)];
sorted_delta_group_rafinc                                           = sortrows(delta_group_rafinc,1,'descend');
sorted_delta_numbers_rafinc                                         = str2double(sorted_delta_group_rafinc(:,1));

% Plot in pathway order
delta_values_rafinc_pway                                            = zeros(length(pway_order_num),1);
for i=1:length(pway_order_num)
     delta_values_rafinc_pway(i)                                    = abs(trapz(allValues_rafinc_mutant_norm(:,pway_order_num(i))) - trapz(allValues_rafinc_wt_norm(:,pway_order_num(i))));
end
delta_group_rafinc_pway                                             = [delta_values_rafinc_pway string(pway_order_num)];

[delta_group_rafinc_pway_nonzero_idx, vals]                         = find(delta_values_rafinc_pway > 0);  


delta_baseline_inc = (delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx)- pway_delta_values(delta_group_pway_nonzero_idx));

% combine plots for base RAF1 vs. 10 fold inc. RAF1
combined_RAF1_settings = [pway_delta_values(delta_group_pway_nonzero_idx) delta_values_rafinc_pway(delta_group_rafinc_pway_nonzero_idx)];
%% 6F
figure
set(gcf,'Position',[100 100 900 300]);
bar(combined_RAF1_settings,1);
hold on
xlim=get(gca,'xlim');
ylim([0 max(combined_RAF1_settings,[],'all')*1.5]);
xlabel('Model Species','FontWeight','bold');
str = {'\int_{0}^{60} |Mutant-WT|(t) dt';'Molecules/cell'};
ylabel(str);
set(gca, 'xTick', 1:length(delta_group_rafinc_pway_nonzero_idx),'xTickLabel',pway_order_group(delta_group_rafinc_pway_nonzero_idx,2),'XTickLabelRotation',45);
title(lowabund_species);
legend('Base RAF1 expression', '[RAF1] * 10')
if savedata ==1
    saveas(gcf,[pwd '/Plots/timeint_diff_baselineRAF1_incRAF1.pdf']);
end
%%
save rasmutant_analyses.mat
