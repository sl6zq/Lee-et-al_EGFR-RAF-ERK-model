%Sloppiness Analysis plots Figs 2, 5, S1 received from co-PI KB 12.4.21

clear all
clc

savedata = false;

paramlist;
%% Load Sloppiness analysis data
twomodesdata        = load('egfr-twomodes.mat');
sloppy_order        = readcell('egfr-parameterorder.txt');
sloppy_cat          = categorical(sloppy_order);
sloppy_cat          = reordercats(sloppy_cat,sloppy_order);
correct_order_mode0 = order_data(twomodesdata.mode_0', sloppy_cat, PLSR_cat_porder);
correct_order_mode1 = order_data(twomodesdata.mode_1', sloppy_cat, PLSR_cat_porder);
%% Plot best-fit values against reference values
bestfit_vals  = params(wanted_param_porder);
% reference values
baselineparam_init_list; %'params' will be defined again
baseline_vals = params(wanted_param_porder);


h1 = figure;
h1.Position = [100 30 500 1100];
for i=1:length(bestfit_vals)
    subplot(10,4,i)
    xline(baseline_vals(i),'--k','LineWidth',2);
    hold on 
    xline(bestfit_vals(i),'--r','LineWidth',2);
    title(names_porder{i}, 'FontSize', 7);
    xlim([min([baseline_vals(i) bestfit_vals(i)])*0.1 max([baseline_vals(i) bestfit_vals(i)])*1.5]);
    %xlim([-0.1 max([baseline_vals(i) bestfit_vals(i)])*1.5]);

    %xlim([min([log(baseline_vals(i)) log(bestfit_vals(i))]) max([log(baseline_vals(i)) log(bestfit_vals(i))])*1.5]);
    check_ph  = ['Best-fit ' names_porder{i} ' = ' num2str(bestfit_vals(i))];
    check_ph2 = ['Baseline ' names_porder{i} ' = ' num2str(baseline_vals(i))];
    %txt = round(bestfit_vals(i)/baseline_vals(i),2);
    fold_diff = bestfit_vals(i)-baseline_vals(i);
    if fold_diff > 1000
        %txt = sprintf('%3.0e',bestfit_vals(i)/baseline_vals(i));
        txt = sprintf('%3.0e', bestfit_vals(i)/ baseline_vals(i));
    else
        txt = round(bestfit_vals(i)/baseline_vals(i),2);
    end
    t = text(max([baseline_vals(i) bestfit_vals(i)])*1.1,1,num2str(txt));
    %t = text(max([log(baseline_vals(i)) log(bestfit_vals(i))])*1.1,1,num2str(txt));
    disp(check_ph)
    disp(check_ph2)

end


h2 = figure;
h2.Position = [100 30 600 1000];
barx = categorical(["Ref.", "Best-fit"]);
barX = reordercats(barx, ["Ref.", "Best-fit"]);
for i=1:length(bestfit_vals)
    subplot(8,5,i)
    bary = [baseline_vals(i) bestfit_vals(i)];
    bar(barX, bary);
    title(names_porder{i}, 'FontSize', 7);
    ylim([min([baseline_vals(i) bestfit_vals(i)])*0.1 max([baseline_vals(i) bestfit_vals(i)])*1.5]);
    %xlim([-0.1 max([baseline_vals(i) bestfit_vals(i)])*1.5]);

    %xlim([min([log(baseline_vals(i)) log(bestfit_vals(i))]) max([log(baseline_vals(i)) log(bestfit_vals(i))])*1.5]);
    check_ph  = ['Best-fit ' names_porder{i} ' = ' num2str(bestfit_vals(i))];
    check_ph2 = ['Baseline ' names_porder{i} ' = ' num2str(baseline_vals(i))];
    %txt = round(bestfit_vals(i)/baseline_vals(i),2);
    fold_diff = bestfit_vals(i)-baseline_vals(i);
    if fold_diff > 1000
        %txt = sprintf('%3.0e',bestfit_vals(i)/baseline_vals(i));
        txt = sprintf('%3.0e', bestfit_vals(i)/ baseline_vals(i));
    else
        txt = round(bestfit_vals(i)/baseline_vals(i),2);
    end
    t = text(max([baseline_vals(i) bestfit_vals(i)])*1.1,1,num2str(txt));
    %t = text(max([log(baseline_vals(i)) log(bestfit_vals(i))])*1.1,1,num2str(txt));
    disp(check_ph)
    disp(check_ph2)

end
if savedata == 1
    saveas(gcf,'sloppycellfit_ref_bargraph.pdf')
end
% log difference
h3 = figure;
h3.Position = [100 30 700 1100];
baseline_vals2 = zeros(length(baseline_vals), 1);
bestfit_vals2  = zeros(length(bestfit_vals), 1);
% Highlight top 5 stiff and sloppy parameters
[sorted_mode0vals, sorted_mode0idx] = sort(abs(correct_order_mode0),'descend');
[sorted_mode1vals, sorted_mode1idx] = sort(abs(correct_order_mode1),'descend');

top_num = 3;
top5stiff_mode0  = PLSR_cat_porder(sorted_mode0idx(1:top_num,1));
top5sloppy_mode0 = PLSR_cat_porder(sorted_mode0idx(end-top_num-1:end,1));
top5stiff_mode1  = PLSR_cat_porder(sorted_mode1idx(1:top_num,1));
top5sloppy_mode1 = PLSR_cat_porder(sorted_mode1idx(end-top_num-1:end,1));
stiff_intersect  = union(top5stiff_mode0, top5stiff_mode1);
sloppy_intersect = union(top5sloppy_mode0, top5sloppy_mode1);
stiff_idx        = find(ismember(PLSR_cat_porder,stiff_intersect));
sloppy_idx       = find(ismember(PLSR_cat_porder,sloppy_intersect));

for i=1:length(bestfit_vals)
    subplot(8,5,i)
    baseline_vals2(i) = log10(baseline_vals(i));
    bestfit_vals2(i)  = log10(bestfit_vals(i));
    xline(baseline_vals2(i),'--k','LineWidth',2);
    hold on 
    xline(bestfit_vals2(i),'--r','LineWidth',2);
    title(names_porder{i}, 'FontSize', 7);
    %xlim([min([baseline_vals2(i) bestfit_vals2(i)]) max([baseline_vals2(i) bestfit_vals2(i)])]);
    xlim([-15 20])
    %xlim([min([log(baseline_vals(i)) log(bestfit_vals(i))]) max([log(baseline_vals(i)) log(bestfit_vals(i))])*1.5]);
    check_ph  = ['Best-fit ' names_porder{i} ' = ' num2str(bestfit_vals2(i))];
    check_ph2 = ['Baseline ' names_porder{i} ' = ' num2str(baseline_vals2(i))];
    %txt = round(bestfit_vals(i)/baseline_vals(i),2);
    fold_diff = abs(bestfit_vals2(i)-baseline_vals2(i));
    if fold_diff > 1000
        %txt = sprintf('%3.0e',bestfit_vals(i)/baseline_vals(i));
        txt = sprintf('%3.0e', bestfit_vals2(i)- baseline_vals2(i));
    else
        %txt = round(bestfit_vals2(i)/baseline_vals2(i),2);
        txt = round(abs(bestfit_vals2(i)-baseline_vals2(i)),3)
    end
    %t = text(max([baseline_vals2(i) bestfit_vals2(i)])*1.3,0.8,num2str(txt));
    %t = text(8,0.8,num2str(txt), 'FontSize',9);
    %hAnnot(i) = annotation('textbox', [i*.2-.2 .6 .1 .1],'String', num2str(txt),'FontSize',9);
    if ismember(i, stiff_idx)
        t = text(20,0.8,'stiff', 'FontSize',9);
    elseif ismember(i, sloppy_idx)
        disp(i)
        t = text(20,0.8,'sloppy', 'FontSize',9);
    end
    disp(check_ph)
    disp(check_ph2)

end

interestingAxes = subplot(8,5,stiff_idx(1));
interestingAxes2 = subplot(8,5,stiff_idx(2));

invisibleAxes = axes('Visible','off','Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'HitTest','off');
invisibleAxes2 = axes('Visible','off','Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'HitTest','off');
border = rectangle('Parent',invisibleAxes,'Position',get(interestingAxes,'Position'),'EdgeColor', [113, 35, 205]/255,...
    'LineWidth',2,'HitTest','off');
%
%border2 = rectangle('Parent',invisibleAxes2,'Position',get(interestingAxes2,'Position'),'EdgeColor', [148, 110, 247]/255,'LineWidth',2,'HitTest','off');
hold off 

if savedata == 1
    saveas(gcf,'sloppycellfit_ref_logscale.pdf')
end


%bar graph of log 10 parameters
h4 = figure;
h4.Position = [100 30 700 1100];
colors_bar  = [0 0 0; 1 1 1];
for i=1:length(bestfit_vals)
    subplot(8,5,i)
    disp(i)
    baseline_vals2(i) = log10(baseline_vals(i));
    bestfit_vals2(i)  = log10(bestfit_vals(i));
    bary = [baseline_vals2(i) bestfit_vals2(i)];
    b = bar(barX, bary);
    for j = 1:numel(b)
        b(j).FaceColor = colors_bar(j,:)
    end
    hold on 
    title(names_porder{i}, 'FontSize', 7);
    ylim([-15 7]);
%     if all([baseline_vals2(i) bestfit_vals2(i)] >= 0)
%         [baseline_vals2(i) bestfit_vals2(i)]
%         ylim([min([baseline_vals2(i) bestfit_vals2(i)])*0.8 max([baseline_vals2(i) bestfit_vals2(i)])*1.2]);
%         
%     num_neg = sum([baseline_vals2(i) bestfit_vals2(i)] <0);
%     elseif num_neg == 1
%         [baseline_vals2(i) bestfit_vals2(i)]
%         ylim([min([baseline_vals2(i) bestfit_vals2(i)])*(2) max([baseline_vals2(i) bestfit_vals2(i)])*10]);
%     elseif all([baseline_vals2(i) bestfit_vals2(i)] <= 0)
%         [baseline_vals2(i) bestfit_vals2(i)]
%         ylim([min([baseline_vals2(i) bestfit_vals2(i)])*(1.5) max([baseline_vals2(i) bestfit_vals2(i)])*(-1.5)]);
%     end
    check_ph  = ['Best-fit ' names_porder{i} ' = ' num2str(bestfit_vals2(i))];
    check_ph2 = ['Baseline ' names_porder{i} ' = ' num2str(baseline_vals2(i))];
    %txt = round(bestfit_vals(i)/baseline_vals(i),2);
    fold_diff = abs(bestfit_vals2(i)-baseline_vals2(i));
    if fold_diff > 1000
        %txt = sprintf('%3.0e',bestfit_vals(i)/baseline_vals(i));
        txt = sprintf('%3.0e', bestfit_vals2(i)- baseline_vals2(i));
    else
        %txt = round(bestfit_vals2(i)/baseline_vals2(i),2);
        txt = round(abs(bestfit_vals2(i)-baseline_vals2(i)),3);
    end

    %hAnnot(i) = annotation('textbox', [i*.4-.2 .6 .1 .1],'String', num2str(txt),'FontSize',9);
    xlims  = 1;
    ylims  = ylim;
    ax = gca;
    axPos  = ax.Position;
    yannot = axPos(2) + axPos(4) * 0.1;
    xannot = axPos(1) + axPos(3) * 0.1;
    %hAnnot(i) = annotation('textbox', [xannot yannot .1 .1],'String', num2str(txt),'FontSize',9, 'Units', 'normalized');
    %hAnnot(i) = text(0.7, 0.1,[num2str(txt)],'FontSize',9, 'Units', 'normalized');
    hAnnot(i) = annotation('textbox', [xannot yannot .1 .1],'String', num2str(txt),'FontSize',9, 'Units', 'normalized');


    if ismember(i, stiff_idx)
        t = text(20,0.8,'stiff', 'FontSize',9);
    elseif ismember(i, sloppy_idx)
        disp(i)
        t = text(20,0.8,'sloppy', 'FontSize',9);
    end
    %disp(check_ph)
    %disp(check_ph2)

end

interestingAxes = subplot(8,5,stiff_idx(1));
interestingAxes2 = subplot(8,5,stiff_idx(2));
interestingAxes3 = subplot(8,5,stiff_idx(3));

invisibleAxes = axes('Visible','off','Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'HitTest','off');
invisibleAxes2 = axes('Visible','off','Position',[0 0 1 1],'XLim',[0 1],'YLim',[0 1],'HitTest','off');

border = rectangle('Parent',invisibleAxes,'Position',get(interestingAxes,'Position'),'EdgeColor', [113, 35, 205]/255,...
    'LineWidth',2,'HitTest','off');
border = rectangle('Parent',invisibleAxes,'Position',get(interestingAxes2,'Position'),'EdgeColor', [113, 35, 205]/255,...
    'LineWidth',2,'HitTest','off');
border = rectangle('Parent',invisibleAxes,'Position',get(interestingAxes3,'Position'),'EdgeColor', [113, 35, 205]/255,...
    'LineWidth',2,'HitTest','off');
hold off 


vals_all = zeros(length(bestfit_vals), 2);
for i=1:length(bestfit_vals)
    vals_all(i,1) = baseline_vals(i);
    vals_all(i,2) = bestfit_vals(i);
end
h2 = figure;
h3.Position = [100 30 700 1000];
barh(vals_all);
legend('Reference', 'Best-fit');
%xl = xline(2328,'-.','Average','DisplayName','Average Sales');

vals_all = log10([max(traj_sormin.RasGTP_best), trapz(traj_sormin.t, traj_sormin.RasGTP_best);...
    max(traj_sorplus.RasGTP_best), trapz(traj_sormin.t, traj_sorplus.RasGTP_best);...
    max(traj_sormin.Raf1_best), trapz(traj_sormin.t, traj_sormin.Raf1_best);...
    max(traj_sorplus.Raf1_best), trapz(traj_sormin.t, traj_sorplus.Raf1_best); ...
    max(traj_sormin.pMEK_best), trapz(traj_sormin.t, traj_sormin.pMEK_best); ...
    max(traj_sormin.pERK_best),  trapz(traj_sormin.t, traj_sormin.pERK_best)]); 

figure
set(gcf,'Position',[100 100 400 350]);
X_cat = categorical({'-sor GTP-bound RAS', '+sor GTP-bound RAS', '-sor membrane RAF1',  '+sor membrane RAF1',  '-sor pMEK', '-sor pERK'});
X_cat = reordercats(X_cat,{'-sor pERK', '-sor pMEK', '+sor membrane RAF1','-sor membrane RAF1','+sor GTP-bound RAS','-sor GTP-bound RAS'});
barh(X_cat, vals_all);
xlim([0 13]);
%set(gca,'yticklabel',["-sor GTP-bound RAS", "+sor GTP-bound RAS", "-sor membrane RAF1",  "+sor membrane RAF1",  "-sor pMEK", "-sor pERK"])
xlabel('log(molecules)'); %set(gca, 'YTick', categorical(vals_all));%ylim([0 max(max(vals_all))*1.2]);
legend('Maximum', 'Average');