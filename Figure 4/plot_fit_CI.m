function plot_fit_CI(all_ci_timeconf, all_ci_output_conf, all_ci_output_mean)

datasets;
xlab        = 'Time (min)';
ylab        = 'Concentration (molec/cell)';
meanlw      = 2;
edgecolor   = 'none';
errorbarcol = 'black';
errorbarlw  = 1;
CI_col_tp   = 0.3;

figure
p1           = fill(all_ci_timeconf,all_ci_output_conf{1}, [1 0.2 0.2],'FaceAlpha', CI_col_tp);
p1.FaceColor = [1, 0.2, 0.2];
p1.EdgeColor = edgecolor;
hold on
p2           = fill(all_ci_timeconf,all_ci_output_conf{2},[0 0.6 0.3],'FaceAlpha', CI_col_tp);
p2.FaceColor = [0 0.6 0.3];
p2.EdgeColor = edgecolor;
plot(min_sor_rastimedata,min_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', [1 0.2 0.2],'LineStyle', edgecolor);
plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor', [0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', edgecolor);
e1           = errorbar(min_sor_rastimedata, min_sor_rasdatanums, min_sor_rasdataerror,'LineStyle', edgecolor);
e1.Color     = errorbarcol;
e1.LineWidth = errorbarlw;
p12          = plot(linspace(0,60),all_ci_output_mean{1},'Color',[1 0.2 0.2],'LineWidth', meanlw);
p22          = plot(linspace(0,60), all_ci_output_mean{2},'Color',[0 0.6 0.3],'LineWidth', meanlw);
xlabel(xlab);
ylabel(ylab);
legend([p12 p22],{'-sor GTP-bound Ras' '+sor GTP-bound Ras'});

figure
p3           = fill(all_ci_timeconf, all_ci_output_conf{3},[1 0.5 0],'FaceAlpha',CI_col_tp);
p3.FaceColor = [1 0.5 0];
p3.EdgeColor = edgecolor;
hold on
p4           = fill(all_ci_timeconf, all_ci_output_conf{4},[0 0.5 1],'FaceAlpha',CI_col_tp);
p4.FaceColor = [0 0.5 1];
p4.EdgeColor = edgecolor;
plot(min_sor_raftimedata,min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.5 0],'LineStyle', edgecolor);
e2           = errorbar(min_sor_raftimedata,min_sor_rafdatanums,min_sor_rafdataerror,'LineStyle','none');
e2.Color     = errorbarcol;
e2.LineWidth = errorbarlw;
plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor', [0 0.5 1], 'MarkerEdgeColor', [0 0.5 1],'LineStyle', edgecolor);
e3           = errorbar(plus_sor_raftimedata,plus_sor_rafdatanums,plus_sor_rafdataerror,'LineStyle','none');
e3.Color     = errorbarcol;
e3.LineWidth = errorbarlw;
p32          = plot(linspace(0,60), all_ci_output_mean{3},'Color',[1 0.5 0],'LineWidth', meanlw);
p42          = plot(linspace(0,60), all_ci_output_mean{4},'Color',[0 0.5 1],'LineWidth', meanlw);
xlabel(xlab);
ylabel(ylab);
legend([p32 p42],{'-sor membrane Raf1' '+sor membrane Raf1'});

figure
p5           = fill(all_ci_timeconf,all_ci_output_conf{5},[0.2 0.8 0.8],'FaceAlpha', CI_col_tp);
p5.FaceColor = [0.2 0.8 0.8];
p5.EdgeColor = edgecolor;
hold on
p6           = fill(all_ci_timeconf,all_ci_output_conf{6},[236/255 0 140/255],'FaceAlpha', CI_col_tp);
p6.FaceColor = [236/255 0 140/255];
p6.EdgeColor = edgecolor;
plot(min_sor_perktimedata,min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor', [0.2 0.8 0.8], 'MarkerEdgeColor', [0.2 0.8 0.8],'LineStyle', edgecolor);
e4           = errorbar(min_sor_perktimedata,min_sor_perkdatanums,min_sor_perkdataerror,'LineStyle', edgecolor);
e4.Color     = errorbarcol;
e4.LineWidth = errorbarlw;
plot(min_sor_pmektimedata,min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor', [236/255 0 140/255], 'MarkerEdgeColor',[236/255 0 140/255],'LineStyle', edgecolor);
e5           = errorbar(min_sor_pmektimedata,min_sor_pmekdatanums,min_sor_pmekdataerror,'LineStyle', edgecolor);
e5.Color     = errorbarcol;
e5.LineWidth = errorbarlw;
p52          = plot(linspace(0,60), all_ci_output_mean{5},'Color',[0.2 0.8 0.8],'LineWidth', meanlw);
p62          = plot(linspace(0,60), all_ci_output_mean{6},'Color',[236/255 0 140/255],'LineWidth', meanlw);
xlabel(xlab);
ylabel(ylab);
legend([p52 p62],{'-sor pERK' '-sor pMEK'});
hold off

end
