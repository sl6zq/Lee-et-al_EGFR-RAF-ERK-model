function [my_remain_params, my_y_max_remain, sorted_aggregate_error] = filter_fnc(number, filter_number, timeSpan, my_other_params, min_sor_total_output, plus_sor_total_output, my_y_max, plotoption)


% Load experimental data
datasample = 50;
datasets;

plus_sor_ras_pred = zeros(length(timeSpan),number);
plus_sor_raf_pred = zeros(length(timeSpan),number);
min_sor_ras_pred = zeros(length(timeSpan),number);
min_sor_raf_pred = zeros(length(timeSpan),number);
min_sor_perk_pred = zeros(length(timeSpan),number);
min_sor_pmek_pred = zeros(length(timeSpan),number);

remain_ind = cell(6,1);
filter_ind = cell(6,1);
array_nums = 1:number;
aggregate_error = zeros(number,1);
aggregate_error_others = zeros(number,5);

for i=1:number
    iter_info = ['Iteration: ' num2str(i)];
    disp(iter_info)
    plus_sor_ras_pred(:,i) = plus_sor_total_output{i}(:,184);
    plus_sor_raf_pred(:,i) = plus_sor_total_output{i}(:,157);
    min_sor_ras_pred(:,i) = min_sor_total_output{i}(:,184);
    min_sor_raf_pred(:,i) = min_sor_total_output{i}(:,157);
    min_sor_perk_pred(:,i) = min_sor_total_output{i}(:,45);
    min_sor_pmek_pred(:,i) = min_sor_total_output{i}(:,21);
    %{
        aggregate_error(i) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_rastimedata,'min_sor_rastimedata') - min_sor_rasdatanums).^2) ./ (length(min_sor_rasdatanums)).^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
            + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_raftimedata,'min_sor_raftimedata') - min_sor_rafdatanums).^2)./ (length(min_sor_rafdatanums)).^2)... %((min_sor_rafdataerror).^2 * length(min_sor_rafdatanums).^2))...
            + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_rastimedata,'plus_sor_rastimedata') - plus_sor_rasdatanums).^2)./ (length(plus_sor_rasdatanums)).^2)...%((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))...
            + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_raftimedata,'plus_sor_raftimedata') - plus_sor_rafdatanums).^2)./ (length(plus_sor_rafdatanums)).^2)...%((plus_sor_rafdataerror).^2 * length(plus_sor_rafdatanums).^2))...
            + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_perktimedata,'min_sor_perktimedata') - min_sor_perkdatanums).^2)./ (length(min_sor_perkdatanums)).^2)...%((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2))...
            + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_pmektimedata,'min_sor_pmektimedata') - min_sor_pmekdatanums).^2./ (length(min_sor_pmekdatanums)).^2)); %((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2)));
    
    %}
    
    %{
    aggregate_error(i) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_rastimedata,'min_sor_rastimedata') - min_sor_rasdatanums).^2) ./ ((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))... %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_raftimedata,'min_sor_raftimedata') - min_sor_rafdatanums).^2)./ ((min_sor_rafdataerror).^2 * length(min_sor_rafdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_rastimedata,'plus_sor_rastimedata') - plus_sor_rasdatanums).^2)./ ((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))... %(length(plus_sor_rasdatanums)).^2)...%((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_raftimedata,'plus_sor_raftimedata') - plus_sor_rafdatanums).^2)./ ((plus_sor_rafdataerror).^2 * length(plus_sor_rafdatanums).^2))...    
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_perktimedata,'min_sor_perktimedata') - min_sor_perkdatanums).^2)./ ((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2))... %(length(min_sor_perkdatanums)).^2)...%((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_pmektimedata,'min_sor_pmektimedata') - min_sor_pmekdatanums).^2./ ((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2))); %(length(min_sor_pmekdatanums)).^2)); %((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2)));
    %}
    
    %{   
    aggregate_error(i) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_rastimedata,'min_sor_rastimedata') - min_sor_rasdatanums).^2) ./ ((min_sor_rasdataerror).^2))... %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_raftimedata,'min_sor_raftimedata') - min_sor_rafdatanums).^2)./ ((min_sor_rafdataerror).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_rastimedata,'plus_sor_rastimedata') - plus_sor_rasdatanums).^2)./ ((plus_sor_rasdataerror).^2))... %(length(plus_sor_rasdatanums)).^2)...%((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_raftimedata,'plus_sor_raftimedata') - plus_sor_rafdatanums).^2)./ ((plus_sor_rafdataerror).^2))...    
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_perktimedata,'min_sor_perktimedata') - min_sor_perkdatanums).^2)./ ((min_sor_perkdataerror).^2))... %(length(min_sor_perkdatanums)).^2)...%((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_pmektimedata,'min_sor_pmektimedata') - min_sor_pmekdatanums).^2./ ((min_sor_pmekdataerror).^2))); %(length(min_sor_pmekdatanums)).^2)); %((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2)));
    %}
    aggregate_error(i) = 1/2 * (sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_rastimedata,'min_sor_rastimedata') - min_sor_rasdatanums).^2) ./ ((min_sor_rasdataerror).^2) ./ (length(min_sor_rasdatanums)))... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_raftimedata,'min_sor_raftimedata') - min_sor_rafdatanums).^2)./ ((min_sor_rafdataerror).^2) ./ (length(min_sor_rafdatanums)))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_rastimedata,'plus_sor_rastimedata') - plus_sor_rasdatanums).^2)./ ((plus_sor_rasdataerror).^2) ./ (length(plus_sor_rasdatanums)))...%((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_raftimedata,'plus_sor_raftimedata') - plus_sor_rafdatanums).^2)./ ((plus_sor_rafdataerror).^2) ./ (length(plus_sor_rafdatanums)))...    
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_perktimedata,'min_sor_perktimedata') - min_sor_perkdatanums).^2)./ ((min_sor_perkdataerror).^2) ./ (length(min_sor_perkdatanums)))...%((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2))...
        + sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_pmektimedata,'min_sor_pmektimedata') - min_sor_pmekdatanums).^2./ ((min_sor_pmekdataerror).^2) ./ (length(min_sor_pmekdatanums))))); %((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2)));
    
    
    
    aggregate_error_others(i,1) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_raftimedata,'min_sor_raftimedata') - min_sor_rafdatanums).^2) ./ ((min_sor_rafdataerror).^2 * length(min_sor_rafdatanums).^2)); %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
    aggregate_error_others(i,2) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_rastimedata,'plus_sor_rastimedata') - plus_sor_rasdatanums).^2)./ ((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2)); %(length(plus_sor_rasdatanums)).^2)...%((plus_sor_rasdataerror).^2 * length(plus_sor_rasdatanums).^2))...
    aggregate_error_others(i,3) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), plus_sor_raftimedata,'plus_sor_raftimedata') - plus_sor_rafdatanums).^2) ./ ((plus_sor_rafdataerror).^2 * length(plus_sor_rafdatanums).^2)); %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
    aggregate_error_others(i,4) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_perktimedata,'min_sor_perktimedata') - min_sor_perkdatanums).^2)./ ((min_sor_perkdataerror).^2 * length(min_sor_perkdatanums).^2)); %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...
    aggregate_error_others(i,5) = sum(((fullEGFR9_onemodel_fit_40(my_other_params(i,:), min_sor_pmektimedata,'min_sor_pmektimedata') - min_sor_pmekdatanums).^2./ ((min_sor_pmekdataerror).^2 * length(min_sor_pmekdatanums).^2))); %(length(min_sor_rasdatanums))^2)... %((min_sor_rasdataerror).^2 * length(min_sor_rasdatanums).^2))...


end
[sorted_aggregate_error, error_index] = sort(aggregate_error,'ascend');
% filter_number = number*1/3;
remain_total = error_index(1:filter_number);
%remain_total = error_index(:,filter_number+1:end);

%{
figure
plot(1:number*1/3, sorted_aggregate_error(:,1:1000));
hold on
xlabel('Parameter set');
ylabel('Aggregate error');
title('Parameter sets with smaller aggregate error');
savefig('aggregate_error_comparison.fig');
%}
plus_sor_ras_remain_final = zeros(length(timeSpan), length(remain_total));
plus_sor_raf_remain_final = zeros(length(timeSpan), length(remain_total));
min_sor_ras_remain_final = zeros(length(timeSpan), length(remain_total));
min_sor_raf_remain_final = zeros(length(timeSpan), length(remain_total));
min_sor_perk_remain_final = zeros(length(timeSpan), length(remain_total));
min_sor_pmek_remain_final = zeros(length(timeSpan), length(remain_total));

for i=1:length(remain_total)
    plus_sor_ras_remain_final(:,i) = plus_sor_total_output{remain_total(i)}(:,184);
    plus_sor_raf_remain_final(:,i) = plus_sor_total_output{remain_total(i)}(:,157);
    min_sor_ras_remain_final(:,i) = min_sor_total_output{remain_total(i)}(:,184);
    min_sor_raf_remain_final(:,i) = min_sor_total_output{remain_total(i)}(:,157);
    min_sor_perk_remain_final(:,i) = min_sor_total_output{remain_total(i)}(:,45);
    min_sor_pmek_remain_final(:,i) = min_sor_total_output{remain_total(i)}(:,21);
end

if strcmp(plotoption, 'on') == 1
    figure
    for i=1:length(remain_total)
        subplot(3,2,1);
        plot(timeSpan,plus_sor_ras_remain_final(:,i))
        hold on
        plot(plus_sor_rastimedata,plus_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(plus_sor_rastimedata, plus_sor_rasdatanums, plus_sor_rasdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('Total GTP-RAS (molec/cell)');
        title('+sor');
        
        subplot(3,2,2);
        plot(timeSpan,plus_sor_raf_remain_final(:,i))
        hold on
        plot(plus_sor_raftimedata,plus_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(plus_sor_raftimedata, plus_sor_rafdatanums, plus_sor_rafdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('Membrane RAF1 (molec/cell)');
        title('+sor');
        
        subplot(3,2,3);
        plot(timeSpan,min_sor_ras_remain_final(:,i))
        hold on
        plot(min_sor_rastimedata, min_sor_rasdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(min_sor_rastimedata, min_sor_rasdatanums, min_sor_rasdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('Total GTP-RAS (molec/cell)');
        title('-sor');
        
        subplot(3,2,4);
        plot(timeSpan,min_sor_raf_remain_final(:,i))
        hold on
        plot(min_sor_raftimedata, min_sor_rafdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(min_sor_raftimedata, min_sor_rafdatanums, min_sor_rafdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('Membrane RAF1 (molec/cell)');
        title('-sor');
        
        
        subplot(3,2,5);
        plot(timeSpan,min_sor_perk_remain_final(:,i))
        hold on
        plot(min_sor_perktimedata, min_sor_perkdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(min_sor_perktimedata, min_sor_perkdatanums, min_sor_perkdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('pERK (molec/cell)');
        title('-sor');
        
        
        subplot(3,2,6);
        plot(timeSpan,min_sor_pmek_remain_final(:,i))
        hold on
        plot(min_sor_pmektimedata, min_sor_pmekdatanums,'Marker', 'o', 'MarkerFaceColor',[0 0.6 0.3], 'MarkerEdgeColor', [0 0.6 0.3],'LineStyle', 'none');
        e1 = errorbar(min_sor_pmektimedata, min_sor_pmekdatanums, min_sor_pmekdataerror,'LineStyle','none');
        e1.Color = 'black';
        e1.LineWidth = 1.4;
        xlabel('time (min)');
        ylabel('pMEK (molec/cell)');
        title('-sor');
        sgtitle('Solutions for only filtered parameter sets (based on all data)');
    end
    savefig('Filtered_slns_alldata_random.fig');
end

%Select model solutions that correspond to the filtered parameter sets
my_remain_params = zeros(length(remain_total),40);
my_y_max_remain = cell(6,1);
for k=1:length(remain_total)
    my_remain_params(k,:) = my_other_params(remain_total(k),:);
    for j=1:6
        my_y_max_remain{j}(k,:) = my_y_max{j}(remain_total(k),:);
    end
end



end