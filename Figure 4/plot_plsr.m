function plot_plsr(PLSR_cat_input, number, sampling, Xloadings, vipScores, PCTVAR, type, title_input)

%Stores necessary PLSR outputs and produce desired plots 
%Inputs: 
    % PLSR_cat_input: names of parameters in the desired order
    % number: number of parameter sets
    % sampling: sampling method used to generate parameter sets
    % Xloadings: predictor loadings 
    % vipScores: VIP scores 
    % PCTVAR: percent variance explained by the model 
    % type: specify the type of desired plot- Xloadings in PC1 or VIP plots
%Outputs: 
    % PLSR parameter loadings in PLS Component 1 
    % or 
    % VIP plots
    
    if strcmp(type, 'Xloadings') == 1
        figure
        bar(PLSR_cat_input,Xloadings(:,1));
        t0 = {['PLSR Principal Component 1 (' num2str(number) ' ' sampling ' parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
        title(t0);
        ylabel('Parameter loadings','FontWeight','bold','FontSize',8);
        set(gcf,'color','w');
        title(title_input);
        ax = gca;
        ax.TitleFontSizeMultiplier = 1;
        set(gca,'FontSize',8);

    end
    if strcmp(type, 'vip') == 1
        [vip_Sorted,vip_index] = sort(vipScores,'descend');
        sorted_vip = vipScores(vip_index)
        high_vip = sorted_vip(sorted_vip > 1);
        high_vip_index = find(sorted_vip > 1);
        
        
        
        f=figure
        f.Position = [100 100 800 200];
        bar(sorted_vip,1)
        ylim([0 max(sorted_vip) + max(sorted_vip)*0.3]);
        hold on
        xlim=get(gca,'xlim');
        plot(xlim,[1 1],'-.k','LineWidth',1);
        title(title_input);
        ylabel('VIP score','FontWeight','bold','FontName','Arial','FontSize',8);
        set(gca, 'xTick', 1:length(PLSR_cat_input),'TickLength', [0 0], 'xTickLabel',PLSR_cat_input(vip_index),'FontName','Arial','FontSize',10);
        xtickangle(90);
        bh = bar(1:numel(sorted_vip),diag(sorted_vip),'stacked','FaceColor', [0/255 155/255 250/255],'BarWidth',1);  
        % Highlight parameters w/ VIP scores > 1
        for k=1:length(high_vip_index)
            bh(k).FaceColor = [227/255 111/255 71/255];
        end

    end    
end
