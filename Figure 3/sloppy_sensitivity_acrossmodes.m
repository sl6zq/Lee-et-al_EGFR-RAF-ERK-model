clear all
clc

figsize = [100 30 600 300];
paramlist;
%% Sensitivity across two modes in Sloppiness analysis
twomodesdata = load('egfr-twomodes.mat');
sloppy_order = readcell('egfr-parameterorder.txt');
sloppy_cat = categorical(sloppy_order);
sloppy_cat = reordercats(sloppy_cat,sloppy_order);
correct_order_mode0 = order_data(twomodesdata.mode_0', sloppy_cat, PLSR_cat_porder);
correct_order_mode1 = order_data(twomodesdata.mode_1', sloppy_cat, PLSR_cat_porder);
correct_order_all   = [correct_order_mode0 correct_order_mode1];

%load eigenvalues of each mode
eigvaldata = load('egfr-logsingval.mat');

h1 = figure
h1.Position = figsize
%bar(PLSR_cat, correct_order_mode0,'k');
bar(PLSR_cat_porder, correct_order_mode0,'k');

h2 = figure;
h2.Position = figsize;
bar(PLSR_cat_porder, correct_order_mode1,'k');
ylabel('Mode 1');

% param sensitivity
param_sens = correct_order_mode0 .^2 * eigvaldata.logs(1) + correct_order_mode1 .^2 * eigvaldata.logs(2);
h3 = figure;
h3.Position = figsize;
bar(PLSR_cat_porder, param_sens,'k');
ylabel('Modes 0 & 1');

%% 6 mdoes 
wanted_sloppy = load('V','V'); %Select the sign of sloppiness analysis eigenvalues (V preferred)
data_sloppy = struct2array(wanted_sloppy);

% SA labels
mylabels;

% Sloppiness labels
labels = readtable('param-order.txt','Format','auto');
labels2 = table2array(labels) 
labels3 = string(labels2); %for reference

% Make sloppycell order follow my order
[data_sloppy_sameorder] = order_data(data_sloppy, labels2, PLSR_cat_porder);

param_sens_all = data_sloppy_sameorder(:,1) .^2 * eigvaldata.logs(1) + data_sloppy_sameorder(:,2) .^2 * eigvaldata.logs(2) + ...
    + data_sloppy_sameorder(:,3) .^2 * eigvaldata.logs(3) + data_sloppy_sameorder(:,4) .^2 * eigvaldata.logs(4) + ...
    + data_sloppy_sameorder(:, 5) .^2 * eigvaldata.logs(5) + data_sloppy_sameorder(:,6) .^2 * eigvaldata.logs(6);

% Compare with PLSR SA VIP 
[params_sens_all_sorted,params_sens_all_sorted_idx] = sort(param_sens_all,'descend');
sorted_param_names     = PLSR_cat_porder(params_sens_all_sorted_idx);
PLSR_cat_sorted_sloppy = categorical(sorted_param_names);
h5 = figure;
h5.Position = figsize;
bar(PLSR_cat_sorted_sloppy, params_sens_all_sorted,'k');
ylabel('Modes 0 ~ 5');

load('diff_RAF1_vip_1xRAF1.mat')
names2 = {'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
names_vec2 = ["kEf" "kcatE" "kpMEK" "knfpSOS" "kBr" "kBf" "kdpERK" "kR1r" "kR1f" "knfpBR1" "kpR1" "kSOSr" "kSOSf" "kRgneslow" "kSon" "kG2SOSr" "kG2SOSf" "kiR1r" "kiR1f" "kSoff" "kpERK" "kRhydro" "Kmgneslow" "kdpMEK" "knfpiR1r" "kdR1r" "kiBr" "kdR1f" "kiBf" "kfpBr" "kdpSOS" "kdp" "kdEr" "kdE,f" "kdpR1" "kG2r" "kG2f" "knfpiBr" "kfpR1r" "kEr"]';
PLSR_cat2 = categorical(names2);
PLSR_cat2 = reordercats(PLSR_cat2,{'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdEf'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'});
[PLSR_vip_sameorder] = order_data(diffraf1_vipScores, PLSR_cat2, PLSR_cat_porder);

h4 = figure;
h4.Position = figsize;
combined_dat = [param_sens_all PLSR_vip_sameorder];
bar(PLSR_cat_porder, combined_dat);
ylabel('Sensitivity');
legend('PLSR SA', 'Sloppiness Analysis');
ylim([0 max(combined_dat, [], 'all')*1.5]);
hold off

h5 = figure;
h5.Position = [100 30 600 80];
h5 = heatmap(param_sens_all', 'CellLabelColor','none','GridVisible','off');
h5.XDisplayLabels = PLSR_cat_porder;
h5.YDisplayLabels = {'Sloppiness Analysis'};
%bwr = @(n)interp1([1 2 3], [90/255 0/255 160/255; 1 1 1; 235/255, 107/255, 35/255], linspace(1, 3, n), 'linear');
bwr = @(n)interp1([1 2 3], [25 23 140; 255 255 255; 201 69 122]./255, linspace(1, 3, n), 'linear');

colormap(bwr(200)); 
set(gca,'FontSize',8);
h5.FontSize = 8;
%saveas(h5, 'Sloppy_sensitivity.pdf');

h6 = figure;
h6.Position = [100 30 600 80];
h6 = heatmap(PLSR_vip_sameorder', 'CellLabelColor','none','GridVisible','off');
h6.XDisplayLabels = PLSR_cat_porder;
h6.YDisplayLabels = {'PLSR SA'};
%bwr = @(n)interp1([1 2 3], [90/255 0/255 160/255; 1 1 1; 235/255, 107/255, 35/255], linspace(1, 3, n), 'linear');
bwr = @(n)interp1([1 2 3], [25 23 140; 255 255 255; 201 69 122]./255, linspace(1, 3, n), 'linear');
colormap(bwr(200)); 
set(gca,'FontSize',8);
h6.FontSize = 8;
%saveas(h6, 'PLSR_SA_VIP.pdf');

%% eigenvalues in percent
raw_vals   = zeros(length(eigvaldata.logs),1);
cumulative = 0;
cumulativesums = zeros(length(eigvaldata.logs),1);
for i=1:length(eigvaldata.logs)
    raw_vals(i) = exp(-i);
    cumulative =  cumulative + raw_vals(i);
    cumulativesums(i) = cumulative;
end
all_modes_eigsum = sum(cumulative);
percent_eigval   = cumulativesums./all_modes_eigsum * 100;


inf_1to6 = zeros(6,1);
for i=1:6 
    inf_1to6(i) = exp(eigvaldata.logs(i));
end
inf_1to6_sum = sum(inf_1to6);

inf_1to2 = zeros(2,1);
for i=1:2 
    inf_1to2(i) = exp(eigvaldata.logs(i));
end
inf_1to2_sum = sum(inf_1to2);

inf_1to41 = zeros(41,1);
for k=1:40 
    inf_1to41(k) = exp(eigvaldata.logs(k));
end
inf_1to41_sum = sum(inf_1to41);

inf_6modes = inf_1to6_sum/inf_1to41_sum * 100
stment = ['6 stiffest modes explain ' num2str(inf_6modes) '% of model behavior'];
disp(stment)
inf_2modes = inf_1to2_sum/inf_1to41_sum * 100
stment2 = ['2 stiffest modes explain ' num2str(inf_2modes) '% of model behavior'];
disp(stment2)



h6 = figure;
h6.Position = [100 30 300 200];
yyaxis left
plot(1:40, eigvaldata.logs,'ko','MarkerSize',4);
ylabel('log(\lambda)');
ylim([min(eigvaldata.logs)*1.1,max(eigvaldata.logs)*1.5]);
yyaxis right
plot(1:40, percent_eigval,'k*','MarkerSize',4);
ylim([0,110]);
ylabel('% fit to data explained');
xlabel('Mode');

h6 = figure;
h6.Position = [100 30 300 200];
yyaxis right
plot(1:40, percent_eigval,'k*','MarkerSize',4);
ylim([0,110]);
yyaxis left
plot(1:40, eigvaldata.logs,'ko','MarkerSize',4);
ylabel('log(\lambda)');
ylabel('% fit to data explained');
