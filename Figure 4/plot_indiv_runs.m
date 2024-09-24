function plot_indiv_runs(all_outputs, all_opt_fits)

% Multiple optimization runs (number of runs = datasample)
datasample = 20;
% Load experimental data (6 sets)
datasets;

min_sor_rasoutput  = all_outputs{1};
plus_sor_rasoutput = all_outputs{2};
min_sor_rafoutput  = all_outputs{3};
plus_sor_rafoutput = all_outputs{4};
min_sor_perkoutput = all_outputs{5};
min_sor_pmekoutput = all_outputs{6};


min_sor_ras_optimalfit  = all_opt_fits{1};
plus_sor_ras_optimalfit = all_opt_fits{2};
min_sor_raf_optimalfit  = all_opt_fits{3};
plus_sor_raf_optimalfit = all_opt_fits{4};
min_sor_perk_optimalfit = all_opt_fits{5};
min_sor_pmek_optimalfit = all_opt_fits{6};

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