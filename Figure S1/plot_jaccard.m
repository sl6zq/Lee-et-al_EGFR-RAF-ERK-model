function [threshold, unique_high_jaccard_all] = plot_jaccard(data_matlab, data_sloppy, names, names_vec)

% Inlcudes pairs in upper right quadrants
ind_matlab = cell(6,6);
ind_sloppy = cell(6,6);
components = 1:6;
unique_high_matlab = cell(6,1);
unique_high_sloppy = cell(6,1);

threshold = 0:0.1:0.6;
for g=1:length(threshold)
    for i = 1:length(components)
        for k = 1:components(i)
            ind_matlab{k,i} = find(abs(data_matlab(:,k)) > threshold(g));
            ind_sloppy{k,i} = find(abs(data_sloppy(:,k)) > threshold(g));
        end
        unique_high_matlab{i} = unique(vertcat(ind_matlab{:,i}));
        unique_high_sloppy{i} = unique(vertcat(ind_sloppy{:,i}));
       
        unique_high_jaccard_all(g,i) = size(intersect(unique_high_matlab{i}, unique_high_sloppy{i}),1) ./ size(union(unique_high_matlab{i}, unique_high_sloppy{i}),1);
    end
end

% figure
% for k=1:length(components)
%     plot(threshold, unique_high_jaccard_all(:,k),'-o');
%     xlabel('Threshold');
%     ylabel('Jaccard index');
%     hold on
%     legend('PC 1', 'PC 1~2', 'PC 1~3', 'PC 1~4', 'PC 1~5', 'PC 1~6','Location','best');
% end
% hold off

for k=1:6
    [values_matlab(:,k), index_matlab(:,k)] = maxk(abs(data_matlab(:,k)),10);
    [values_sloppy(:,k), index_sloppy(:,k)] = maxk(abs(data_sloppy(:,k)),10);
    top_data_matlab{k} = [names_vec(index_matlab(:,k)) values_matlab(:,k)];
    top_data_sloppy{k} = [names_vec(index_sloppy(:,k)) values_sloppy(:,k)];
    jaccard{k} = size(intersect(top_data_matlab{k}(:,1), top_data_sloppy{k}(:,1)),1) ./ size(union(top_data_matlab{k}(:,1), top_data_sloppy{k}(:,1)),1);
end
% Inlcudes all pairs
figure
for k=1:6
    a(k) = subplot(2,3,k);
    scatter(abs(data_matlab(:,k)),abs(data_sloppy(:,k)),'filled','MarkerFaceColor', [63/255, 97/255, 215/255], 'MarkerEdgeColor', [63/255, 97/255, 215/255]);
    labelpoints(abs(data_matlab(:,k)),abs(data_sloppy(:,k)), names, 'SE', 0.2, 1);
    hold on
    %h = refline;
    xlabel('MATLAB');
    ylabel('SloppyCell');
    str2 = {['J = ', num2str(round(jaccard{1,k},3))]};
    %xlim2 = xlim; ylim2 = ylim;
    xlim2 = 0.8; ylim2 = 0.8;
    xlim([0 xlim2]);
    ylim([0 ylim2]);
    %text(xlim2 - xlim2/4, ylim2 - ylim2/7.5, str2);
    %b(k) = annotation('textbox', 'String', str2, 'HorizontalAlignment', 'Right', 'FitBoxToText','on');
    %text(0.2, 0.2, str2);
    title1 = ['PC', num2str(k), ' / Mode ', num2str(k-1)];
    title(title1);
    sgtitle('All parameters(Jaccard Index for top 10 parameters)');
end

end
