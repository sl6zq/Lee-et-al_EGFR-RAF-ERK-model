function plot_correlation(number, sampling, data_matlab, data_sloppy, names)
    
    fig = figure;
    fig.Position = [100 100 620 220];
    coeff = zeros(size(data_matlab,2),1);
    p_val = zeros(size(data_matlab,2),1);
    val_threshold = 0.5; 
    for i=1:size(data_matlab,2)
        s1 = subplot(2,3,i);
        %plot(data_matlab(:,i),data_sloppy(:,i),'o','MarkerFaceColor', [62/255, 154/255, 77/255],'MarkerEdgeColor',[62/255, 154/255, 77/255]);
        plot(data_matlab(:,i),data_sloppy(:,i),'o','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',4);
        hold on
        l = lsline ;
        set(l,'LineWidth', 1, 'Color', 'r', 'LineStyle',"--");
        line([0 0.5], [0.5 0.5], 'Color', 'k', 'LineStyle',":");
        line([0.5 0.5], [0 0.5], 'Color', 'k', 'LineStyle',":");

        mat_df = data_matlab(:,i);
        sloppy_df = data_sloppy(:,i);
        matlab_strong = find(abs(data_matlab(:,i)) > val_threshold);
        sloppy_strong = find(abs(data_sloppy(:,i)) > val_threshold);
        both_strong = union(matlab_strong, sloppy_strong);
        diffmat = data_matlab(:,i) - data_sloppy(:,i);
        %T = table(names, data_matlab(:,i), data_sloppy(:,i), diffmat(:,i));
        %diffmat = sortrows(T,'Var4');
        %topfive = diffmat(1:5)
        %[idx_topfive, x] = ismember(topfive, diffmat)
        labelpoints(data_matlab(both_strong,i)+0.03,data_sloppy(both_strong,i)+0.05, names(both_strong), 'SE', 0.2, 1, 'FontSize', 7); %south east
        [coeff(i), p_val(i)] = corr(data_matlab(:,i), data_sloppy(:,i));
        S = sprintf('p = %0.2e', p_val(i));
        S = regexprep(S, 'e\+?(-?\d+)', ' x 10^{$1}');

        str1 = {['R = ', num2str(round(corr(data_matlab(:,i), data_sloppy(:,i)),2))]};
        str2 = {S};
        s1.XLim = [0 1];
        s1.YLim = [0 1];
        text(s1.XLim(:,2)-0.52,s1.YLim(:,2)-0.1, str1, 'FontSize',7);
        text(s1.XLim(:,2)-0.52,s1.YLim(:,2)-0.2, str2, 'FontSize',7);
        title_str = ['PC ', num2str(i), '/Mode ', num2str(i-1)];
        title(title_str,'FontWeight','bold');

        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabeltext = {['Sloppiness Analysis eigenvalues']};
        xlabeltext = {['PLSR SA Loadings']};
        ylabel(han,ylabeltext);
        xlabel(han,xlabeltext);
        filename = ['pearsons ', num2str(number), sampling '.fig'];
        %saveas(fig, filename);
        hold off
    end
end