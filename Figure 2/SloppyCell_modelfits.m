%Model fits based on SloppyCell-fitted parameter values, received from co-PI KB 12.4.21

clear all
clc

paramlist;
datasets;

%% A: Model fits 
traj_sormin  = load('ens-trajectories-sorminus');
traj_sorplus = load('ens-trajectories-sorplus');

target      = max(traj_sormin.Raf1_best) * 0.05; %time point where 5% of peak membrane RAF1
closest_idx = find(traj_sormin.Raf1_best == interp1(traj_sormin.Raf1_best,traj_sormin.Raf1_best,target,'nearest')); 
new_colors  = {[0 0.5 1],[1 0.5 0],[0 0.5 1],[1 0.5 0],[0 0.5 1],[0 0.5 1]};

figure
set(gcf,'Position',[100 100 400 320]);
m1 = plot(traj_sormin.t,traj_sormin.RasGTP_mu,'--','color', new_colors{1}, 'LineWidth', 2);
hold on
plot(traj_sormin.t,traj_sormin.RasGTP_best,'-','color', new_colors{1}, 'LineWidth', 2);
p1 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.RasGTP_mu-traj_sormin.RasGTP_sigma  fliplr(traj_sormin.RasGTP_mu + traj_sormin.RasGTP_sigma)], new_colors{1}, 'FaceAlpha',0.3);
p1.EdgeColor = 'none';
m2 = plot(traj_sorplus.t,traj_sorplus.RasGTP_mu,'--','color', new_colors{2}, 'LineWidth', 2);
p2 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [traj_sorplus.RasGTP_mu-traj_sorplus.RasGTP_sigma  fliplr(traj_sorplus.RasGTP_mu + traj_sorplus.RasGTP_sigma)], new_colors{2}, 'FaceAlpha',0.3);
plot(traj_sorplus.t,traj_sorplus.RasGTP_best,'-','color', new_colors{2}, 'LineWidth', 2);
p2.EdgeColor = 'none';
xlabel('Time (min)');
ylabel('Concentration (molec/cell)');
d1 = plot(min_sor_rastimedata,min_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','LineStyle', 'none');
d2 = plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k','LineStyle', 'none');
e1 = errorbar(min_sor_rastimedata,min_sor_rasdatanums,min_sor_rasdataerror,'LineStyle','none');
e1.Color = 'black';
e1.LineWidth = 1;
d3 = xline(traj_sormin.t(closest_idx),'--k'); 
legend([d1 d2 d3],{'-sor GTP-bound Ras' '+sor GTP-bound Ras' '5% of active RAF1 remains'},'FontSize',9);
ylim([0 max(traj_sormin.RasGTP_best) * 1.2]);
xlim([0 62]);
ax = gca;
hold off

figure
set(gcf,'Position',[100 100 400 320]);
m3 = plot(traj_sormin.t,traj_sormin.Raf1_mu,'--','color', new_colors{3}, 'LineWidth', 2);
hold on
plot(traj_sormin.t,traj_sormin.Raf1_best,'-','color', new_colors{3}, 'LineWidth', 2);
p3 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.Raf1_mu-traj_sormin.Raf1_sigma  fliplr(traj_sormin.Raf1_mu + traj_sormin.Raf1_sigma)], new_colors{3}, 'FaceAlpha',0.3);
p3.EdgeColor = 'none';
m4 = plot(traj_sormin.t,traj_sorplus.Raf1_mu,'--','color', new_colors{4}, 'LineWidth', 2);
plot(traj_sormin.t,traj_sorplus.Raf1_best,'-','color', new_colors{4}, 'LineWidth', 2);
p4 = patch([traj_sorplus.t fliplr(traj_sorplus.t)], [traj_sorplus.Raf1_mu-traj_sorplus.Raf1_sigma  fliplr(traj_sorplus.Raf1_mu + traj_sorplus.Raf1_sigma)], new_colors{4}, 'FaceAlpha',0.3);
p4.EdgeColor = 'none';
d4 = plot(min_sor_raftimedata,min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', new_colors{3}, 'MarkerEdgeColor', new_colors{3},'LineStyle', 'none');
e2 = errorbar(min_sor_raftimedata,min_sor_rafdatanums,min_sor_rafdataerror,'LineStyle','none');
e2.Color = new_colors{3};
e2.LineWidth = 1;
d5 = plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', new_colors{4}, 'MarkerEdgeColor', new_colors{4},'LineStyle', 'none');
e3 = errorbar(plus_sor_raftimedata,plus_sor_rafdatanums,plus_sor_rafdataerror,'LineStyle','none');
e3.Color = new_colors{4};
e3.LineWidth = 1;
xlabel('Time (min)');
ylabel('Concentration (molec/cell)');
ylim([0 Inf]);
d6 = xline(traj_sormin.t(closest_idx),'--k'); 
legend([d4 d5 d6],{'-sor membrane Raf1' '+sor membrane Raf1' '5% of active RAF1 remains'});
xlim([0 62]);ylim([-500 max(plus_sor_rafdatanums) * 1.2]);
hold off

figure
set(gcf,'Position',[100 100 400 320]);
m5 = plot(traj_sormin.t,traj_sormin.pMEK_mu,'--','color', new_colors{5}, 'LineWidth', 2);
hold on
plot(traj_sormin.t,traj_sormin.pMEK_best,'-','color', new_colors{5}, 'LineWidth', 2);
p5 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.pMEK_mu-traj_sormin.pMEK_sigma  fliplr(traj_sormin.pMEK_mu + traj_sormin.pMEK_sigma)], new_colors{5}, 'FaceAlpha',0.3);
p5.EdgeColor = 'none';
d7 = plot(min_sor_pmektimedata,min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor',new_colors{5}, 'MarkerEdgeColor',new_colors{5},'LineStyle', 'none');
e5 = errorbar(min_sor_pmektimedata,min_sor_pmekdatanums,min_sor_pmekdataerror,'LineStyle','none');
e5.Color = new_colors{5};
e5.LineWidth = 1;
xlabel('Time (min)');
ylabel('Concentration (molec/cell)');
ylim([0 max(min_sor_perkdatanums) * 1.2]);
xlim([0 62]);
d8 = xline(traj_sormin.t(closest_idx),'--k'); 
legend([d7 d8],{'-sor pMEK' '5% of active RAF1 remains'});
hold off

figure
set(gcf,'Position',[100 100 400 320]);
plot(traj_sormin.t,traj_sormin.pERK_best,'-','color', new_colors{6}, 'LineWidth', 2);
hold on
m6 = plot(traj_sormin.t,traj_sormin.pERK_mu,'--','color', new_colors{6}, 'LineWidth', 2);
p6 = patch([traj_sormin.t fliplr(traj_sormin.t)], [traj_sormin.pERK_mu-traj_sormin.pERK_sigma  fliplr(traj_sormin.pERK_mu + traj_sormin.pERK_sigma)], new_colors{6}, 'FaceAlpha',0.3);
p6.EdgeColor = 'none';
d9 = plot(min_sor_perktimedata,min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor', new_colors{6}, 'MarkerEdgeColor',new_colors{6},'LineStyle', 'none');
e4 = errorbar(min_sor_perktimedata,min_sor_perkdatanums,min_sor_perkdataerror,'LineStyle','none');
e4.Color = new_colors{6};
e4.LineWidth = 1;
d10 = xline(traj_sormin.t(closest_idx),'--k'); 
legend([d9 d10],{'-sor pERK' '5% of active RAF1 remains'});
ylim([0 max(min_sor_perkdatanums) * 1.2]);ax = gca;
xlim([0 62]);

%% Compare membrane RAF1 vs. membrane SOS trajectories
[T1, ~, ~, params_minsor, ~, allValues_rand]   = fullEGFR9_onemodel(0:1:60, yinit, params','min_sor','');
% membrane SOS
%figure
%plot(T1, allValues_rand(:,157)) %Raf1_pm_param
%hold on
%plot(T1, allValues_rand(:,44)) %EG2SOS
%legend('Membrane RAF1', 'Membrane SOS');
%hold off

percent_SOS_memb = ['% SOS membrane-bound: ' num2str(max(allValues_rand(:,44))/yinit(34) * 100)];
disp(percent_SOS_memb)
%% Average % deviation of confidence interval bounds and exp. error
min_sor_pmek_int_lb  = interp1(traj_sormin.t,(traj_sormin.pMEK_mu-traj_sormin.pMEK_sigma),min_sor_pmektimedata);
min_sor_pmek_int_ub  = interp1(traj_sormin.t,(traj_sormin.pMEK_mu+traj_sormin.pMEK_sigma),min_sor_pmektimedata);

meandata_all  = {min_sor_rasdatanums; plus_sor_rasdatanums; min_sor_rafdatanums; plus_sor_rafdatanums; min_sor_pmekdatanums; min_sor_perkdatanums};
errordata_all = {min_sor_rasdataerror; plus_sor_rasdataerror; min_sor_rafdataerror; plus_sor_rafdataerror; min_sor_pmekdataerror; min_sor_perkdataerror};
meanpred_all  = {traj_sormin.RasGTP_mu; traj_sorplus.RasGTP_mu; traj_sormin.Raf1_mu; traj_sorplus.Raf1_mu; traj_sormin.pMEK_mu; traj_sormin.pERK_mu};
meantime_all  = {traj_sormin.t; traj_sorplus.t; traj_sormin.t; traj_sorplus.t; traj_sormin.t; traj_sormin.t};
meansigma_all = {traj_sormin.RasGTP_sigma; traj_sorplus.RasGTP_sigma; traj_sormin.Raf1_sigma; traj_sorplus.Raf1_sigma; traj_sormin.pMEK_sigma; traj_sormin.pERK_sigma};
timedat_all   = {min_sor_rastimedata; plus_sor_rastimedata; min_sor_raftimedata; plus_sor_raftimedata; min_sor_pmektimedata; min_sor_perktimedata};

pctdiff       = cell(length(meandata_all),1);
for k=1:length(meandata_all)
    int_lb    = interp1(meantime_all{k},(meanpred_all{k} - meansigma_all{k}),timedat_all{k});
    int_ub    = interp1(meantime_all{k},(meanpred_all{k} + meansigma_all{k}),timedat_all{k});
    %figure
    %plot(meantime_all{k},meanpred_all{k}-meansigma_all{k},'--',meantime_all{k},meanpred_all{k}+meansigma_all{k},'--',timedat_all{k},int_ub,':.');
    %hold on
    %e4 = errorbar(timedat_all{k},meandata_all{k},errordata_all{k},'LineStyle','none');
    % how far off is upper bound of predicted values from data error
    int_data_ub     = meandata_all{k} + errordata_all{k};
    int_data_lb     = meandata_all{k} - errordata_all{k};
    int_ub_pctdiff  = (int_ub - int_data_ub)./int_data_ub * 100;
    int_lb_pctdiff  = (int_lb - int_data_lb)./int_data_lb * 100;
    comb_pctdiff    = [int_ub_pctdiff int_lb_pctdiff];
    pctdiff_lb_ub   = mean(comb_pctdiff,2);
    pctdiff{k}      = mean(pctdiff_lb_ub);
end

all_pctdiff    = vertcat(pctdiff{:});
avg_all_pctdff = mean(all_pctdiff,1);

% average RMSE 
bestfit_rmse        = cell(length(meandata_all),1);
for k=1 %:length(meandata_all)
    int_mean        = interp1(meantime_all{k},(meanpred_all{k}),timedat_all{k});
    %figure
    %plot(timedat_all{k},int_mean,'-k');
    %hold on
    %e4 = errorbar(timedat_all{k},meandata_all{k},errordata_all{k},'LineStyle','none');
    %plot(timedat_all{k},meandata_all{k},'o');
    % how far off is best-fit values from mean data (RMSE)
    MSE             = (int_mean - meandata_all{k}).^2;
    bestfit_rmse{k} = sqrt(sum(MSE)./length(meandata_all{k}));
end
all_rmse  = vertcat(bestfit_rmse{:});
avg_rmse  = mean(all_rmse,1);

%% RMSE (of mean fits) and standard dev
% extract data time points
[t_val, t_id] = find(traj_sormin.t==min_sor_rastimedata);
% solve ODEs again for specific data time points
rmsd_all = (traj_sormin.RasGTP_mu - min_sor_rasdatanums);


%% B: Max. & time-averaged model outputs (Fig. 3B) for -sor, + sor species
paramlist;

[T1, ~, ~, params_minsor, ~, allValues_minsor]   = fullEGFR9_onemodel(traj_sormin.t, yinit, params,'min_sor', ''); 
[~, ~, ~, params_plussor, ~, allValues_plussor]  = fullEGFR9_onemodel(traj_sormin.t, yinit, params,'plus_sor', ''); 

traj_sormin_pEGFR  = allValues_minsor(:,41)  + allValues_minsor(:,25)  + allValues_minsor(:,44); %'E', 'EG2', 'EG2SOS'
traj_plussor_pEGFR = allValues_plussor(:,41) + allValues_plussor(:,25) + allValues_plussor(:,44); %'E' ,' EG2', 'EG2SOS'

max_vals  = ([max(traj_sormin_pEGFR), max(traj_plussor_pEGFR), ...
    max(traj_sormin.RasGTP_best), max(traj_sorplus.RasGTP_best), ...
    max(traj_sormin.Raf1_best), max(traj_sorplus.Raf1_best), ...
    max(traj_sormin.pMEK_best), max(traj_sormin.pERK_best)]);
tint_vals = ([trapz(traj_sormin.t, traj_sormin_pEGFR), trapz(traj_sormin.t, traj_sormin_pEGFR)/max(traj_sormin.t), ...
    trapz(traj_sormin.t, traj_sormin.RasGTP_best), trapz(traj_sormin.t, traj_sorplus.RasGTP_best)/max(traj_sormin.t), ...
    trapz(traj_sormin.t, traj_sormin.Raf1_best), trapz(traj_sormin.t, traj_sorplus.Raf1_best)/max(traj_sormin.t), ...
    trapz(traj_sormin.t, traj_sormin.pMEK_best), trapz(traj_sormin.t, traj_sormin.pERK_best)/max(traj_sormin.t)]);
max_t = max(traj_sormin.t);
vals_all = ([max(traj_sormin_pEGFR), trapz(traj_sormin.t, traj_sormin_pEGFR)/max_t;...
    max(traj_plussor_pEGFR), trapz(traj_sormin.t, traj_plussor_pEGFR)/max_t;...
    max(traj_sormin.RasGTP_best), trapz(traj_sormin.t, traj_sormin.RasGTP_best)/max_t;...
    max(traj_sorplus.RasGTP_best), trapz(traj_sormin.t, traj_sorplus.RasGTP_best)/max_t;...
    max(traj_sormin.Raf1_best), trapz(traj_sormin.t, traj_sormin.Raf1_best)/max_t;...
    max(traj_sorplus.Raf1_best), trapz(traj_sormin.t, traj_sorplus.Raf1_best)/max_t; ...
    max(traj_sormin.pMEK_best), trapz(traj_sormin.t, traj_sormin.pMEK_best)/max_t; ...
    max(traj_sormin.pERK_best),  trapz(traj_sormin.t, traj_sormin.pERK_best)/max_t]); 

figure
set(gcf,'Position',[100 100 400 350]);
X_cat = categorical({'-sor pEGFR', '+sor pEGFR', '-sor GTP-bound RAS', '+sor GTP-bound RAS', '-sor membrane RAF1',  '+sor membrane RAF1',  '-sor pMEK', '-sor pERK'});
X_cat = reordercats(X_cat,{'-sor pERK', '-sor pMEK', '+sor membrane RAF1','-sor membrane RAF1','+sor GTP-bound RAS','-sor GTP-bound RAS','+sor pEGFR', '-sor pEGFR'});
barh(X_cat, vals_all);
xlabel('molecules');
legend('Maximum', 'Average');
