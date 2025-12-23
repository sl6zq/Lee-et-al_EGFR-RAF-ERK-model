clc
clear all

numsim = 198;

% Read files file1.txt through file20.txt. Files are in the current directory.

for k = 1 : numsim
	textFileName = ['mathmodel_',num2str(k-1), '.txt'];
	if isfile(textFileName)
		fid = fopen(textFileName, 'rt');
		tline = fgetl(fid);
        headers = strsplit(tline, ' ');     %a cell array of strings
        %48 species for stochastic
        datacell{k} = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
        fclose(fid);
        timescan = cell2mat(datacell{1}(1));    %as a numeric array
	else
		fprintf('File %s does not exist.\n', textFileName);
	end
end

figure
for k = 1 : numsim
    Ras_GTP(:,k) = (datacell{k}{:,20} + (2.0 * (datacell{k}{:,27} + datacell{k}{:,28} + datacell{k}{:,26} + datacell{k}{:,23} + datacell{k}{:,30} + datacell{k}{:,29} + datacell{k}{:,34} + datacell{k}{:,38} + datacell{k}{:,33} + datacell{k}{:,41} + datacell{k}{:,31} + datacell{k}{:,32})) + datacell{k}{:,24} + datacell{k}{:,22} + datacell{k}{:,25} + datacell{k}{:,21} + datacell{k}{:,36} + datacell{k}{:,40} + datacell{k}{:,35} + datacell{k}{:,39} + datacell{k}{:,37});
    Raf1_pm(:,k) = (datacell{k}{:,25} + (2.0 * datacell{k}{:,23}) + datacell{k}{:,26} + datacell{k}{:,22} + datacell{k}{:,21} + (2.0 * datacell{k}{:,27}) + datacell{k}{:,28} + datacell{k}{:,30} + (2.0 * datacell{k}{:,29}) + datacell{k}{:,40} + datacell{k}{:,36} + datacell{k}{:,41} + datacell{k}{:,31} + datacell{k}{:,32});
    pMEK(:,k) = cell2mat(datacell{k}(:,8));
    pERK(:,k) = cell2mat(datacell{k}(:,11));
    
    subplot(2,2,1)
    plot(timescan,Ras_GTP(:,k));
    title('Total GTP-bound Ras');
    hold on 
        
    subplot(2,2,2)
    plot(timescan,Raf1_pm(:,k));
    title('Membrane-bound Raf');
    hold on 
        
    subplot(2,2,3)
    plot(timescan,pMEK(:,k));
    title('pMEK');
    hold on 
        
    subplot(2,2,4)
    plot(timescan,pERK(:,k));
    title('pERK');
    sgtitle('-sorafenib');
    hold on 
end
hold off

%deterministic
for k = 1 : numsim
	textFileName_det = [num2str(k-1), '.txt'];
	if isfile(textFileName_det)
		fid = fopen(textFileName_det, 'rt');
		tline = fgetl(fid);
        headers = strsplit(tline, ' ');     %a cell array of strings
        %47 species for stochastic
        datacell_det{k} = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
        fclose(fid);
        timescan = cell2mat(datacell_det{1}(1));    %as a numeric array
	else
		fprintf('File %s does not exist.\n', textFileName_det);
	end
end

my_y_Ras_GTP = zeros(numsim,1);
my_y_Raf1_pm = zeros(numsim,1);
my_y_pMEK = zeros(numsim,1);
my_y_pERK = zeros(numsim,1);

figure
for k = 1 : numsim
    Ras_GTP_det(:,k) = (datacell_det{k}{:,19} + (2.0 * (datacell_det{k}{:,26} + datacell_det{k}{:,27} + datacell_det{k}{:,25} + datacell_det{k}{:,22} + datacell_det{k}{:,29} + datacell_det{k}{:,28} + datacell_det{k}{:,33} + datacell_det{k}{:,37} + datacell_det{k}{:,32} + datacell_det{k}{:,40} + datacell_det{k}{:,30} + datacell_det{k}{:,31})) + datacell_det{k}{:,23} + datacell_det{k}{:,21} + datacell_det{k}{:,24} + datacell_det{k}{:,20} + datacell_det{k}{:,35} + datacell_det{k}{:,39} + datacell_det{k}{:,34} + datacell_det{k}{:,38} + datacell_det{k}{:,36});
    Raf1_pm_det(:,k) = (datacell_det{k}{:,24} + (2.0 * datacell_det{k}{:,22}) + datacell_det{k}{:,25} + datacell_det{k}{:,21} + datacell_det{k}{:,20} + (2.0 * datacell_det{k}{:,26}) + datacell_det{k}{:,27} + datacell_det{k}{:,29} + (2.0 * datacell_det{k}{:,28}) + datacell_det{k}{:,39} + datacell_det{k}{:,35} + datacell_det{k}{:,40} + datacell_det{k}{:,30} + datacell_det{k}{:,31});
    pMEK_det(:,k) = cell2mat(datacell_det{k}(:,7));
    pERK_det(:,k) = cell2mat(datacell_det{k}(:,10));
    
    Ras_GTP_residual(:,k) = Ras_GTP_det(:,k) - Ras_GTP(:,k);
    Raf1_pm_residual(:,k) = Raf1_pm_det(:,k) - Raf1_pm(:,k);
    pMEK_residual(:,k) = pMEK_det(:,k) - pMEK(:,k);
    pERK_residual(:,k) = pERK_det(:,k) - pERK(:,k);
    
    
    Ras_GTP_var_residual(:,k) = var(Ras_GTP_residual(:,k))';
    Raf1_pm_var_residual(:,k) = var(Raf1_pm_residual(:,k))';
    pMEK_var_residual(:,k) = var(pMEK_residual(:,k))';
    pERK_var_residual(:,k) = var(pERK_residual(:,k))';
    
    new_max_Ras_GTP = max(Ras_GTP_det(:,k));
    my_y_Ras_GTP(k) = new_max_Ras_GTP;
    new_max_Raf1_pm = max(Raf1_pm_det(:,k));
    my_y_Raf1_pm(k) = new_max_Raf1_pm;
    new_max_pMEK = max(pMEK_det(:,k));
    my_y_pMEK(k) = new_max_pMEK;
    new_max_pERK = max(pERK_det(:,k));
    my_y_pERK(k) = new_max_pERK;
    
    
    subplot(2,2,1)
    plot(timescan,Ras_GTP_det(:,k));
    title('Total GTP-bound Ras');
    hold on 
        
    subplot(2,2,2)
    plot(timescan,Raf1_pm_det(:,k));
    title('Membrane-bound Raf');
    hold on 
        
    subplot(2,2,3)
    plot(timescan,pMEK_det(:,k));
    title('pMEK');
    hold on 
        
    subplot(2,2,4)
    plot(timescan,pERK_det(:,k));
    title('pERK');
    sgtitle('-sorafenib');
    hold on 
end
hold off

figure
subplot(2,2,1)
plot(1:numsim,Ras_GTP_var_residual)
hold on
title('GTP-bound Ras');
xlabel('Parameter array');
ylabel('Variance of residual');
%set stochasticity threshold to 0.5*10^5 
%plot(xlim,[0.1*10^5  0.1*10^5 ],'-.')
hold off
subplot(2,2,2)
plot(1:numsim,Raf1_pm_var_residual)
hold on
title('membrane-bound Raf');
xlabel('Parameter array');
ylabel('Variance of residual');
%set stochasticity threshold to 0.5*10^5 
%plot(xlim,[0.1*10^5  0.1*10^5 ],'-.')
hold off
subplot(2,2,3)
plot(1:numsim,pMEK_var_residual)
hold on
title('pMEK');
xlabel('Parameter array');
ylabel('Variance of residual');
%set stochasticity threshold to 0.5*10^5 
%plot(xlim,[0.1*10^5  0.1*10^5 ],'-.')
hold off
subplot(2,2,4)
plot(1:numsim,pERK_var_residual)
hold on
title('pERK');
xlabel('Parameter array');
ylabel('Variance of residual');
%set stochasticity threshold to 0.5*10^5 
%plot(xlim,[0.1*10^5  0.1*10^5 ],'-.')
sgtitle('Degree of stochasticity by each parameter array');
hold off


%Get parameter array text files

%PLSR with output as max
paramdata=readtable('param_array2.dat','Delimiter',',','ReadVariableNames',false)
paramdata = table2cell(paramdata)
%b=regexp(paramdata,'\d+(\.)?(\d+)?','match')
%b=regexp(paramdata,'\d+(\.)?(\d+)?','match')
b = regexp(paramdata, '(?<name>\w+)=(?<value>[^,]+)', 'names');
all_values_as_numeric = zeros(size(b,1),size(b,2));
for m=1:size(b,1)
    for n=1:size(b,2)
        all_values_as_numeric(m,n) = str2double( {b{m,n}.value} );
    end
end

number = numsim;
param_size = size(all_values_as_numeric);
observations = param_size(1);
my_log_params = zeros(param_size);
for o = 1:observations
    for p = 1:size(all_values_as_numeric,2)
        if all_values_as_numeric(o,p) == 0
            my_log_params(o,p) = all_values_as_numeric(o,p);
        else
            my_log_params(o,p) = log10(all_values_as_numeric(o,p));
        end
    end
end
plsr_z_score_x = zeros(number, size(all_values_as_numeric,2)); %-size_r when removing aka one iteration
for p = 1:size(all_values_as_numeric,2)
    x = my_log_params(:,p);
    Zsc = @(x) (x - mean(x))./std(x);   % Z-score function
    Zx = Zsc(x);
    plsr_z_score_x(:,p) = Zx;
end

 
%{
plsr_z_score_y = zeros(number, size(var_residual,2)); %-size_r when removing aka one iteration
for p = 1:size(var_residual,2)
    y = var_residual(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_y(:,p) = Zy;
end
%}
plsr_z_score_Ras_GTP = zeros(number, size(my_y_Ras_GTP,2)); %-size_r when removing aka one iteration
for p = 1:size(my_y_Ras_GTP,2)
    y = my_y_Ras_GTP(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_Ras_GTP(:,p) = Zy;
end
plsr_z_score_Raf1_pm = zeros(number, size(my_y_Raf1_pm,2)); %-size_r when removing aka one iteration
for p = 1:size(my_y_Raf1_pm,2)
    y = my_y_Raf1_pm(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_Raf1_pm(:,p) = Zy;
end
plsr_z_score_pMEK = zeros(number, size(my_y_pMEK,2)); %-size_r when removing aka one iteration
for p = 1:size(my_y_pMEK,2)
    y = my_y_pMEK(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_pMEK(:,p) = Zy;
end
plsr_z_score_pERK = zeros(number, size(my_y_pERK,2)); %-size_r when removing aka one iteration
for p = 1:size(my_y_pERK,2)
    y = my_y_pERK(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_pERK(:,p) = Zy;
end

%plsr_z_score_y = [plsr_z_score_Ras_GTP plsr_z_score_Raf1_pm plsr_z_score_pMEK plsr_z_score_pERK];
plsr_z_score_y = plsr_z_score_Ras_GTP;

ncomp = 10;
[n,p] = size(plsr_z_score_y);
C1 = cvpartition(numsim, 'LeaveOut');

[Xloadings,YL,XS,YS,BETA,PCTVAR,PLSmsep,stats] = plsregress(plsr_z_score_x, plsr_z_score_y, ncomp,'CV',C1);
Q2Y = 1- PLSmsep(2,2:end)/sum(sum((plsr_z_score_y-mean(plsr_z_score_y)).^2)./size(plsr_z_score_y,1));
R2Y = cumsum(PCTVAR(2,1:end));
W0 = bsxfun(@rdivide,stats.W,sqrt(sum(stats.W.^2,1)));
% Calculate the product of summed squares of XS and YL
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
% Calculate VIP scores for NCOMP components
vipScores = sqrt(size(Xloadings,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));

names = {'k_{B,f}'; 'k_{B,r}'; 'k_{catE}'; 'k_{dE,f}'; 'k_{dE,r}'; 'k_{dp}'; 'k_{dpERK}'; 'k_{dpMEK}'; 'k_{dpR1}'; 'k_{dpSOS}'; 'k_{dR1,f}'; 'k_{dR1,r}'; 'k_{E,f}'; 'k_{E,r}'; 'k_{fpBr}'; 'k_{fpR1,r}'; 'k_{G2,f}'; 'k_{G2,r}'; 'k_{G2SOS,f}'; 'k_{G2SOS,r}'; 'k_{iB,f}'; 'k_{iB,r}';
'k_{iR1,f}'; 'k_{iR1,r}'; 'K_{mgneslow}'; 'k_{nfpBR1}'; 'k_{nfpiB,r}'; 'k_{nfpiR1,r}'; 'k_{nfpSOS}'; 'k_{pERK}'; 'k_{pMEK}';
'k_{pR1}'; 'k_{R1,f}'; 'k_{R1,r}'; 'k_{Rgneslow}'; 'k_{Rhydro}'; 'k_{Soff}'; 'k_{Son}'; 'k_{SOS,f}'; 'k_{SOS,r}'}

PLSR_cat = categorical(names);

PLSR_cat = reordercats(PLSR_cat,{'k_{B,f}'; 'k_{B,r}'; 'k_{catE}'; 'k_{dE,f}'; 'k_{dE,r}'; 'k_{dp}'; 'k_{dpERK}'; 'k_{dpMEK}'; 'k_{dpR1}'; 'k_{dpSOS}'; 'k_{dR1,f}'; 'k_{dR1,r}'; 'k_{E,f}'; 'k_{E,r}'; 'k_{fpBr}'; 'k_{fpR1,r}'; 'k_{G2,f}'; 'k_{G2,r}'; 'k_{G2SOS,f}'; 'k_{G2SOS,r}'; 'k_{iB,f}'; 'k_{iB,r}';
'k_{iR1,f}'; 'k_{iR1,r}'; 'K_{mgneslow}'; 'k_{nfpBR1}'; 'k_{nfpiB,r}'; 'k_{nfpiR1,r}'; 'k_{nfpSOS}'; 'k_{pERK}'; 'k_{pMEK}';
'k_{pR1}'; 'k_{R1,f}'; 'k_{R1,r}'; 'k_{Rgneslow}'; 'k_{Rhydro}'; 'k_{Soff}'; 'k_{Son}'; 'k_{SOS,f}'; 'k_{SOS,r}'});




figure
bar(PLSR_cat,Xloadings(:,1));
t1 = {['PLSR Principal Component 1 for Max. Outputs (' num2str(numsim) ' log-uniform random parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
%title(t1);
title('PLSR Principal Component 1');
ylabel('Parameter loadings','FontWeight','bold');
set(gcf,'color','w');
ax = gca;
ax.XAxis.FontSize = 9.5;
ax.YAxis.FontSize = 12;
set(gcf, 'Position',  [10, 10, 700, 500]);
str1 = {['R^2Y: ',num2str(round(R2Y(:,1)),3)],['Q^2Y: ',num2str(round(Q2Y(:,1),3))]};
annotation('textbox', [0.77, 0.8, 0.1, 0.1], 'String', str1, 'FitBoxToText','on','HorizontalAlignment','center');


figure
bar(PLSR_cat,vipScores);
%t1 = {['PLSR Principal Component 1 for Max. Outputs (' num2str(numsim) ' log-uniform random parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
%title(t1);
title('Maximum GTP-bound Ras');
ylabel('VIP Score','FontWeight','bold');
set(gcf,'color','w');
ax = gca;
ax.XAxis.FontSize = 11.2;
ax.YAxis.FontSize = 12;
set(gcf, 'Position',  [10, 10, 800, 500]);
str1 = {['R^2Y: ',num2str(round(R2Y(:,1)),3)],['Q^2Y: ',num2str(round(Q2Y(:,1),3))]};
hold on
plot(xlim,[1 1],'-.')

figure
plot(1:numsim,Ras_GTP_var_residual,'LineWidth',1.2)
hold on
title('GTP-bound Ras');
xlabel('Parameter array','FontSize',12);
ylabel('Variance of residual','FontSize',12);
%set stochasticity threshold to 0.5*10^5 
plot(xlim,[0.1*10^5  0.1*10^5 ],'-.','LineWidth',1.2)
hold off

high_stoch_params_index  = find(Ras_GTP_var_residual > 0.2*10^5);
high_stoch_params_values = Ras_GTP_var_residual(high_stoch_params_index);
remain_stoch_params_index = find(Ras_GTP_var_residual < 0.2*10^5);

all_values_as_numeric_remain_stoch_params = all_values_as_numeric;
all_values_as_numeric_remain_stoch_params(high_stoch_params_index,:) = [];  % remove highly stochastic parameter arrays

param_size_remain = size(all_values_as_numeric_remain_stoch_params);
observations = param_size_remain(1);
my_log_remain_stoch_params = zeros(param_size_remain);
for o = 1:observations
    for p = 1:size(all_values_as_numeric_remain_stoch_params,2)
        if all_values_as_numeric_remain_stoch_params(o,p) == 0
            my_log_remain_stoch_params(o,p) = all_values_as_numeric_remain_stoch_params(o,p);
        else
            my_log_remain_stoch_params(o,p) = log10(all_values_as_numeric_remain_stoch_params(o,p));
        end
    end
end
plsr_z_score_remain_stoch_params = zeros(param_size_remain(1), size(all_values_as_numeric_remain_stoch_params,2)); %-size_r when removing aka one iteration
for p = 1:size(all_values_as_numeric_remain_stoch_params,2)
    x_remain_stoch_params = my_log_remain_stoch_params(:,p);
    Zsc_remain_stoch_params = @(x_remain_stoch_params) (x_remain_stoch_params - mean(x_remain_stoch_params))./std(x_remain_stoch_params);   % Z-score function
    Zx_remain_stoch_params = Zsc(x_remain_stoch_params);
    plsr_z_score_remain_stoch_params(:,p) = Zx_remain_stoch_params;
end

plsr_z_score_Ras_GTP_remain = plsr_z_score_Ras_GTP;
plsr_z_score_Ras_GTP_remain(high_stoch_params_index,:) = []; 

C2 = cvpartition(param_size_remain(1), 'LeaveOut');
[Xloadings_Ras_GTP_remain,YL_Ras_GTP_remain,XS_Ras_GTP_remain,YS_Ras_GTP_remain,BETA_Ras_GTP_remain,PCTVAR_Ras_GTP_remain,PLSmsep_Ras_GTP_remain,stats_Ras_GTP_remain] = plsregress(plsr_z_score_remain_stoch_params, plsr_z_score_Ras_GTP_remain, ncomp,'CV',C2);
Q2Y_Ras_GTP_remain = 1- PLSmsep_Ras_GTP_remain(2,2:end)/sum(sum((plsr_z_score_Ras_GTP_remain-mean(plsr_z_score_Ras_GTP_remain)).^2)./size(plsr_z_score_Ras_GTP_remain,1));
R2Y_Ras_GTP_remain = cumsum(PCTVAR_Ras_GTP_remain(2,1:end));
W0_Ras_GTP_remain = bsxfun(@rdivide,stats.W,sqrt(sum(stats_Ras_GTP_remain.W.^2,1)));
% Calculate the product of summed squares of XS and YL
sumSq_Ras_GTP_remain = sum(XS.^2,1).*sum(YL_Ras_GTP_remain.^2,1);
% Calculate VIP scores for NCOMP components
vipScores_Ras_GTP_remain = sqrt(size(Xloadings_Ras_GTP_remain,1) * sum(bsxfun(@times,sumSq_Ras_GTP_remain,W0_Ras_GTP_remain.^2),2) ./ sum(sumSq_Ras_GTP_remain,2));

figure
plot(1:ncomp,cumsum(100*PCTVAR(2,:)),'-bo');
hold on
plot(1:ncomp,cumsum(100*PCTVAR_Ras_GTP_remain(2,:)),'-ro');
xlabel('Number of PLS components','FontSize',12);
ylabel('Percent Variance Explained in Y','FontSize',12);
tplsr = [num2str(number) ' parameter arrays'];
t2plsr = [num2str(param_size_remain(1)) ' parameter arrays'];
legend(tplsr, t2plsr,'Location','Best','FontSize',12);

%Use of KS statistics (test for 198 parameter arrays from
%fullEGFR9_MathModel_stochastic_for_batch_params_SL_copy)

%read parameter array 
for k = 1 : numsim
	textFileName_param_sens = ['param_sens',num2str(k-1), '.txt'];
	if isfile(textFileName_param_sens)
		fid_param_sens = fopen(textFileName_param_sens, 'rt');
		tline = fgetl(fid_param_sens);
        headers = strsplit(tline, ' ');     %a cell array of strings
        %47 species for stochastic
        datacell_det{k} = textscan(fid_param_sens, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
        fclose(fid_param_sens);
        timescan = cell2mat(datacell_det{1}(1));    %as a numeric array
	else
		fprintf('File %s does not exist.\n', textFileName_param_sens);
	end
end

[h0,p0,k2stat0] = kstest(x_default, x_0, pd_default, pd_0)
[h1,p1,k2stat1] = kstest(x_default, x_1, pd_default, pd_1)
[h2,p2,k2stat2] = kstest(x_default, x_2, pd_default, pd_2)
[h3,p3,k2stat3] = kstest(x_default, x_3, pd_default, pd_3)
[h4,p4,k2stat4] = kstest(x_default, x_4, pd_default, pd_4)
[h5,p5,k2stat5] = kstest(x_default, x_5, pd_default, pd_5)
[h6,p6,k2stat6] = kstest(x_default, x_6, pd_default, pd_6)
[h7,p7,k2stat7] = kstest(x_default, x_7, pd_default, pd_7)
[h8,p8,k2stat8] = kstest(x_default, x_8, pd_default, pd_8)
[h9,p9,k2stat9] = kstest(x_default, x_9, pd_default, pd_9)

k2stat_total= [k2stat0 k2stat1 k2stat2 k2stat3 k2stat4 k2stat5 k2stat6 k2stat7 k2stat8 k2stat9];
h_total= [h0 h1 h2 h3 h4 h5 h6 h7 h8 h9]

figure
b = bar(1:10,k2stat_total);
xlabel('Parameter Array');
ylabel('KSP');
title('SOS');
xtips = find(h_total==0);
ytips = b(1).YEndPoints(:,xtips);
%labels1 = string(b(1:10).YData);
labels1 = ["H=0", "H=0", "H=0", "H=0", "H=0"];
text(xtips,ytips,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

figure
bar(1:10,h_total,'r');
xlabel('Parameter Array');
ylabel('KS test decision');
title('SOS');

figure
for i=1:length(high_stoch_params_index)
    plot(timescan,Ras_GTP(:,high_stoch_params_index(i)));
    hold on
    xlabel('time(min)');
    ylabel('GTP-bound Ras');
    title('Model evaluation with high stochasticity (> 0.2*10^5) parameter arrays');
end

figure
for i=1:length(remain_stoch_params_index)
    plot(timescan,Ras_GTP(:,remain_stoch_params_index(i)));
    hold on
    xlabel('time(min)');
    ylabel('GTP-bound Ras');
    title('Model evaluation with low stochasticity (< 0.2*10^5) parameter arrays');
end

[max_stoch_params_value,max_stoch_params_index]  = max(pERK_var_residual);
[min_stoch_params_value,min_stoch_params_index]  = min(pERK_var_residual);

figure
plot(timescan,pERK(:,max_stoch_params_index),'r','LineWidth',1.3);
hold on
plot(timescan,pERK_det(:,max_stoch_params_index),'--r','LineWidth',1.3);
plot(timescan,pERK(:,min_stoch_params_index),'b','LineWidth',1.3);
plot(timescan,pERK_det(:,min_stoch_params_index),'--b','LineWidth',1.3);
xlabel('time(min)');
ylabel('pERK (molec/cell)');
legend('Stochastic w/ max. stochasticity param array', 'Deterministic w/ max. stochasticity param array', 'Stochastic w/ min. stochasticity param array ','Deterministic w/ min. stochasticity param array ');




%Use of KS statistics (test for 10 parameter arrays from
%fullEGFR13_paramarraytest_newparam for ERK_count)


[h0,p0,k2stat0] = kstest(x_default, x_0, pd_default, pd_0)
[h1,p1,k2stat1] = kstest(x_default, x_1, pd_default, pd_1)
[h2,p2,k2stat2] = kstest(x_default, x_2, pd_default, pd_2)
[h3,p3,k2stat3] = kstest(x_default, x_3, pd_default, pd_3)
[h4,p4,k2stat4] = kstest(x_default, x_4, pd_default, pd_4)
[h5,p5,k2stat5] = kstest(x_default, x_5, pd_default, pd_5)
[h6,p6,k2stat6] = kstest(x_default, x_6, pd_default, pd_6)
[h7,p7,k2stat7] = kstest(x_default, x_7, pd_default, pd_7)
[h8,p8,k2stat8] = kstest(x_default, x_8, pd_default, pd_8)
[h9,p9,k2stat9] = kstest(x_default, x_9, pd_default, pd_9)

k2stat_total= [k2stat0 k2stat1 k2stat2 k2stat3 k2stat4 k2stat5 k2stat6 k2stat7 k2stat8 k2stat9]

figure
bar(1:10,k2stat_total);
xlabel('Parameter Array');
ylabel('KSP');
title('SOS');


residual0 = pMEK0_det{1} - pMEK0{1};
residual1 = pMEK1_det{1} - pMEK1{1};
residual2 = pMEK2_det{1} - pMEK2{1};
residual3 = pMEK3_det{1} - pMEK3{1};
residual4 = pMEK4_det{1} - pMEK4{1};
residual5 = pMEK5_det{1} - pMEK5{1};
residual6 = pMEK6_det{1} - pMEK6{1};
residual7 = pMEK7_det{1} - pMEK7{1};
residual8 = pMEK8_det{1} - pMEK8{1};
residual9 = pMEK9_det{1} - pMEK9{1};

totalresidual = [residual0 residual1 residual2 residual3 residual4 residual5 residual6 residual7 residual8 residual9];
var_residual = var(totalresidual)';

%datacell = cell(numsim,49);
for i = 1:2
    fid = fopen(['mathmodel_',sprintf('%03d',num2str(i)),'.txt'], 'rt');
    tline = fgetl(fid);
    headers = strsplit(tline, ' ');     %a cell array of strings
    datacell = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
    fclose(fid);
    timescan = datacell{1};    %as a numeric array
end


names = {'k_{E,f}';
'k_{catE}';
'k_{pMEK}';
'k_{nfpSOS}';
'k_{Br}';
'k_{Bf}';
'k_{dpERK}';
'k_{R1r}';
'k_{R1f}';
'k_{nfpBR1}';
'k_{on}';
'k_{pR1}';
'k_{SOSr}';
'k_{SOSf}';
'k_{Rgneslow}';
'k_{Son}';
'k_{G2SOSr}';
'k_{G2SOSf}';
'k_{iR1r}';
'k_{iR1f}';
'k_{Soff}';
'k_{p}';
'k_{pERK}';
'k_{Rhydro}';
'K_{mgneslow}';
'k_{dpMEK}';
'k_{nfpiR1r}';
'k_{dR1r}';
'k_{iBr}';
'k_{dR1f}';
'k_{iBf}';
'k_{fpBr}';
'k_{off}';
'k_{dpSOS}';
'k_{dp}';
'k_{dEr}';
'k_{dEf}';
'k_{dpR1}';
'k_{G2r}';
'k_{G2f}';
'k_{nfpiBr}';
'k_{fpR1r}';
'k_{Er'};


PLSR_cat = categorical(names);
PLSR_cat = reordercats(PLSR_cat,{'k_{E,f}';
'k_{catE}';
'k_{pMEK}';
'k_{nfpSOS}';
'k_{Br}';
'k_{Bf}';
'k_{dpERK}';
'k_{R1r}';
'k_{R1f}';
'k_{nfpBR1}';
'k_{on}';
'k_{pR1}';
'k_{SOSr}';
'k_{SOSf}';
'k_{Rgneslow}';
'k_{Son}';
'k_{G2SOSr}';
'k_{G2SOSf}';
'k_{iR1r}';
'k_{iR1f}';
'k_{Soff}';
'k_{p}';
'k_{pERK}';
'k_{Rhydro}';
'K_{mgneslow}';
'k_{dpMEK}';
'k_{nfpiR1r}';
'k_{dR1r}';
'k_{iBr}';
'k_{dR1f}';
'k_{iBf}';
'k_{fpBr}';
'k_{off}';
'k_{dpSOS}';
'k_{dp}';
'k_{dEr}';
'k_{dEf}';
'k_{dpR1}';
'k_{G2r}';
'k_{G2f}';
'k_{nfpiBr}';
'k_{fpR1r}';
'k_{Er'});

%{
fid0 = fopen('mathmodel_000.txt', 'rt');
tline0 = fgetl(fid0);
datacell0 = textscan(fid0, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK0 = datacell0(:,8);
pMEK0{1}(end,:) = [];
timescan = datacell0{1}; 
timescan(end,:) = [];
fclose(fid0);

fid1 = fopen('001.txt', 'rt');
tline1 = fgetl(fid1);
datacell1 = textscan(fid1, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK1 = datacell1(:,8);
pMEK1{1}(end,:) = [];
fclose(fid1);

fid2 = fopen('002.txt', 'rt');
tline2 = fgetl(fid2);
datacell2 = textscan(fid2, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK2 = datacell2(:,8);
pMEK2{1}(end,:) = [];
fclose(fid2);

fid3 = fopen('003.txt', 'rt');
tline3 = fgetl(fid3);
datacell3 = textscan(fid3, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK3 = datacell3(:,8);
pMEK3{1}(end,:) = [];
fclose(fid3);

fid4 = fopen('004.txt', 'rt');
tline4 = fgetl(fid4);
datacell4 = textscan(fid4, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK4 = datacell4(:,8);
pMEK4{1}(end,:) = [];
fclose(fid4);

fid5 = fopen('005.txt', 'rt');
tline5 = fgetl(fid5);
datacell5 = textscan(fid5, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK5 = datacell5(:,8);
pMEK5{1}(end,:) = [];
fclose(fid5);

fid6 = fopen('006.txt', 'rt');
tline6 = fgetl(fid6);
datacell6 = textscan(fid6, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK6 = datacell6(:,8);
pMEK6{1}(end,:) = [];
fclose(fid6);

fid7 = fopen('007.txt', 'rt');
tline7 = fgetl(fid7);
datacell7 = textscan(fid7, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK7 = datacell7(:,8);
pMEK7{1}(end,:) = [];
fclose(fid7);

fid8 = fopen('008.txt', 'rt');
tline8 = fgetl(fid8);
datacell8 = textscan(fid8, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK8 = datacell8(:,8);
pMEK8{1}(end,:) = [];
fclose(fid8);

fid9 = fopen('009.txt', 'rt');
tline9 = fgetl(fid9);
datacell9 = textscan(fid9, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK9 = datacell9(:,8);
pMEK9{1}(end,:) = [];
fclose(fid9);
%}
pMEKtotal = cell2mat([pMEK0 pMEK1 pMEK2 pMEK3 pMEK4 pMEK5 pMEK6 pMEK7 pMEK8 pMEK9]); 

figure
for k=1:10
    plot(timescan,pMEKtotal(:,k));
    xlabel('time(min)');
    ylabel('pMEK (molecules/cell)');
    title('Stochastic');
    hold on
end
hold off


fid0_det = fopen('000_det.txt', 'rt');
tline0_det = fgetl(fid0_det);
datacell0_det = textscan(fid0_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK0_det = datacell0_det(:,8);
timescan0_det = datacell0_det{1}; 
fclose(fid0_det);

fid1_det = fopen('001_det.txt', 'rt');
tline1_det = fgetl(fid1_det);
datacell1_det = textscan(fid1_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK1_det = datacell1_det(:,8);
timescan1_det = datacell1_det{1}; 
fclose(fid1_det);

fid2_det = fopen('002_det.txt', 'rt');
tline2_det = fgetl(fid2_det);
datacell2_det = textscan(fid2_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK2_det = datacell2_det(:,8);
timescan2_det = datacell2_det{1}; 
fclose(fid2_det);

fid3_det = fopen('003_det.txt', 'rt');
tline3_det = fgetl(fid3_det);
datacell3_det = textscan(fid3_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK3_det = datacell3_det(:,8);
timescan3_det = datacell3_det{1}; 
fclose(fid3_det);

fid4_det = fopen('004_det.txt', 'rt');
tline4_det = fgetl(fid4_det);
datacell4_det = textscan(fid4_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK4_det = datacell4_det(:,8);
timescan4_det = datacell4_det{1}; 
fclose(fid4_det);

fid5_det = fopen('005_det.txt', 'rt');
tline5_det = fgetl(fid5_det);
datacell5_det = textscan(fid5_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK5_det = datacell5_det(:,8);
timescan5_det = datacell5_det{1}; 
fclose(fid5_det);

fid6_det = fopen('006_det.txt', 'rt');
tline6_det = fgetl(fid6_det);
datacell6_det = textscan(fid6_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK6_det = datacell6_det(:,8);
timescan6_det = datacell6_det{1}; 
fclose(fid6_det);

fid7_det = fopen('007_det.txt', 'rt');
tline7_det = fgetl(fid7_det);
datacell7_det = textscan(fid7_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK7_det = datacell7_det(:,8);
timescan7_det = datacell7_det{1}; 
fclose(fid7_det);

fid8_det = fopen('008_det.txt', 'rt');
tline8_det = fgetl(fid8_det);
datacell8_det = textscan(fid8_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
pMEK8_det = datacell8_det(:,8);
timescan8_det = datacell8_det{1}; 
fclose(fid8_det);

fid9_det = fopen('009_det.txt', 'rt');
tline9_det = fgetl(fid9_det);
datacell9_det = textscan(fid9_det, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',' '); %split data by species (specis per cell)
timescan9_det = datacell9_det{1}; 
pMEK9_det = datacell9_det(:,8);
fclose(fid9_det);

%pMEK_det_total = cell2mat([pMEK0_det pMEK1_det pMEK2_det pMEK3_det pMEK4_det pMEK5_det pMEK6_det pMEK7_det pMEK8_det pMEK9_det]); 

figure
plot(timescan0_det,pMEK0_det{1});
hold on
plot(timescan1_det,pMEK1_det{1});
plot(timescan2_det,pMEK2_det{1});
plot(timescan3_det,pMEK3_det{1});
plot(timescan4_det,pMEK4_det{1});
plot(timescan5_det,pMEK5_det{1});
plot(timescan6_det,pMEK6_det{1});
plot(timescan7_det,pMEK7_det{1});
plot(timescan8_det,pMEK8_det{1});
plot(timescan9_det,pMEK9_det{1});
xlabel('time(min)');
ylabel('pMEK (molecules/cell)');
title('Deterministic');
hold off

figure
plot(timescan0_det,pMEK0_det{1},'r','LineWidth', 1.2);
hold on
plot(timescan,pMEK0{1},'r');
legend('Deterministic','Stochastic');
xlabel('time(min)');
ylabel('pMEK (molecules/cell)');
hold off



residual0 = pMEK0_det{1} - pMEK0{1};
residual1 = pMEK1_det{1} - pMEK1{1};
residual2 = pMEK2_det{1} - pMEK2{1};
residual3 = pMEK3_det{1} - pMEK3{1};
residual4 = pMEK4_det{1} - pMEK4{1};
residual5 = pMEK5_det{1} - pMEK5{1};
residual6 = pMEK6_det{1} - pMEK6{1};
residual7 = pMEK7_det{1} - pMEK7{1};
residual8 = pMEK8_det{1} - pMEK8{1};
residual9 = pMEK9_det{1} - pMEK9{1};

totalresidual = [residual0 residual1 residual2 residual3 residual4 residual5 residual6 residual7 residual8 residual9];
var_residual = var(totalresidual)';

figure
stem(timescan0_det,residual0);
xlabel('time(min)');
ylabel('Residuals');
hold off

figure
stem(timescan5_det,residual5);
xlabel('time(min)');
ylabel('Residuals');
hold off

paramdata=readtable('fullEGFR12param_array.txt','Delimiter',',','ReadVariableNames',false)
paramdata = table2cell(paramdata)
%b=regexp(paramdata,'\d+(\.)?(\d+)?','match')
%b=regexp(paramdata,'\d+(\.)?(\d+)?','match')
b = regexp(paramdata, '(?<name>\w+)=(?<value>[^,]+)', 'names');

all_values_as_numeric = zeros(size(b,1),size(b,2));
for m=1:size(b,1)
    for n=1:size(b,2)
        all_values_as_numeric(m,n) = str2double( {b{m,n}.value} );
    end
end
%out=str2double([b{:}])

%{
fid_param = fopen('fullEGFR12param_array.txt', 'rt');
datacell_param = textscan(fid_param, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter',','); %split data by species (specis per cell)
%raw_varied_params = ; 
%}
number = 10;
param_size = size(all_values_as_numeric);
observations = param_size(1);
my_log_params = zeros(param_size);
for o = 1:observations
    for p = 1:size(all_values_as_numeric,2)
        if all_values_as_numeric(o,p) == 0
            my_log_params(o,p) = all_values_as_numeric(o,p);
        else
            my_log_params(o,p) = log10(all_values_as_numeric(o,p));
        end
    end
end
plsr_z_score_x = zeros(number, size(all_values_as_numeric,2)); %-size_r when removing aka one iteration
for p = 1:size(all_values_as_numeric,2)
    x = my_log_params(:,p);
    Zsc = @(x) (x - mean(x))./std(x);   % Z-score function
    Zx = Zsc(x);
    plsr_z_score_x(:,p) = Zx;
end

plsr_z_score_y = zeros(number, size(var_residual,2)); %-size_r when removing aka one iteration
for p = 1:size(var_residual,2)
    y = var_residual(:,p);
    Zsc = @(y) (y - mean(y))./std(y);   % Z-score function
    Zy = Zsc(y)  ;
    plsr_z_score_y(:,p) = Zy;
end

ncomp = 5;
[n,p] = size(plsr_z_score_y);
C1 = cvpartition(number, 'LeaveOut');

[Xloadings,YL,XS,YS,BETA,PCTVAR,PLSmsep,stats] = plsregress(plsr_z_score_x, plsr_z_score_y, ncomp,'CV',C1);
Q2Y = 1- PLSmsep(2,2:end)/sum(sum((plsr_z_score_y-mean(plsr_z_score_y)).^2)./size(plsr_z_score_y,1));
R2Y = cumsum(PCTVAR(2,1:end));
W0 = bsxfun(@rdivide,stats.W,sqrt(sum(stats.W.^2,1)));
% Calculate the product of summed squares of XS and YL
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
% Calculate VIP scores for NCOMP components
vipScores = sqrt(size(Xloadings,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));



figure
bar(PLSR_cat,Xloadings(:,1));
t1 = {['PLSR Principal Component 1 for Max. Outputs (' num2str(number) ' log-uniform random parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
%title(t1);
title('PLSR Principal Component 1');
ylabel('Parameter loadings','FontWeight','bold');
set(gcf,'color','w');
ax = gca;
ax.XAxis.FontSize = 9.5;
ax.YAxis.FontSize = 12;
set(gcf, 'Position',  [10, 10, 700, 500]);
str1 = {['R^2Y: ',num2str(round(R2Y(:,1)),3)],['Q^2Y: ',num2str(round(Q2Y(:,1),3))]};
annotation('textbox', [0.77, 0.8, 0.1, 0.1], 'String', str1, 'FitBoxToText','on','HorizontalAlignment','center');


figure
bar(PLSR_cat,vipScores);
t1 = {['PLSR Principal Component 1 for Max. Outputs (' num2str(number) ' log-uniform random parameter sets)'];[num2str(PCTVAR(2,1)*100) '% variance in Y ,' num2str(PCTVAR(1,1)*100) '% variance in X explained']};
%title(t1);
%title('PLSR Principal Component 1');
ylabel('VIP Score','FontWeight','bold');
set(gcf,'color','w');
ax = gca;
ax.XAxis.FontSize = 11.2;
ax.YAxis.FontSize = 12;
set(gcf, 'Position',  [10, 10, 800, 500]);
str1 = {['R^2Y: ',num2str(round(R2Y(:,1)),3)],['Q^2Y: ',num2str(round(Q2Y(:,1),3))]};
hold on
plot(xlim,[1 1],'-.')
%annotation('textbox', [0.77, 0.8, 0.1, 0.1], 'String', str1, 'FitBoxToText','on','HorizontalAlignment','center');

%conduct KS test of cdf of the distributions at last time point
%[h,p, ks2stat] = kstest2(cdf1, cdf2)
