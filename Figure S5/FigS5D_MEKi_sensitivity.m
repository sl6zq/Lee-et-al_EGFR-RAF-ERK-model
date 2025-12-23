clear all 
clc
% Multiple optimization runs (number of runs = datasample)
datasample = 20;
% Load experimental data (6 sets)
datasets;

%% Define parameters
rng(1) %random seed for SA
paramlist; % load best-fit values from SloppyCell
%baselineparam_init_list; %baseline parameters 
sampling  = 'random'; %sampling method for SA

% ODE param
num_Param = 3000; % at least 2 for cv partition, num_Param >= 122
timeSpan  = 0:1:60;
% Maximum number of PLS components in PLSR model
ncomp     = 40;
%% Test sensitivity of pERK to MEKi
% Are low RAF-expressing cells more sensitive to MEKi?
% Low RAF
k_MEKi    = linspace(1, 0, 20);
%k_MEKi    = linspace(100, 0, 20);
%k_MEKi    = [0.001 0.01 0.1 1 10 100];
%k_MEKi    = [0.0001 0.001 0.01 0.1 1];
%k_MEKi    = linspace(1, 100, 20);


pERK_cell  = cell(length(k_MEKi),1);
MEK_cell   = cell(length(k_MEKi),1);
rasmut_pERK_cell  = cell(length(k_MEKi),1);
rasmut_MEK_cell   = cell(length(k_MEKi),1);
max_pERK   = zeros(length(k_MEKi),1);
max_MEK    = zeros(length(k_MEKi),1);
rasmut_max_pERK   = zeros(length(k_MEKi),1);
rasmut_max_MEK    = zeros(length(k_MEKi),1);
tint_pERK         = zeros(length(k_MEKi),1);
tint_MEK          = zeros(length(k_MEKi),1);
rasmut_tint_pERK  = zeros(length(k_MEKi),1);
rasmut_tint_MEK   = zeros(length(k_MEKi),1);

paramlist_krasmutant; %16-fold decrease in kRhydro

for i=1:length(k_MEKi)
    [T,~,~,params_minsor,allNames,allValues_wt]             = fullEGFR9_onemodel_MEKi(timeSpan, yinit, params, 'min_sor', 'no', k_MEKi(i));
    [~,~,~,~,~,allValues_rasmut]                            = fullEGFR9_onemodel_MEKi(timeSpan, yinit, params_mutant, 'min_sor', 'no', k_MEKi(i));
    pERK_cell{i} = allValues_wt(:,45); 
    MEK_cell{i} = allValues_wt(:,28); 
    
    rasmut_pERK_cell{i} = allValues_wt(:,45); 
    rasmut_MEK_cell{i} = allValues_wt(:,28); 
    
    max_pERK(i) = max(pERK_cell{i});
    max_MEK(i)  = max(MEK_cell{i});
    
    rasmut_max_pERK(i) = max(rasmut_pERK_cell{i});
    rasmut_max_MEK(i)  = max(rasmut_MEK_cell{i});
    
    tint_pERK(i) = trapz(pERK_cell{i});
    tint_MEK(i)  = trapz(MEK_cell{i});
    
    rasmut_tint_pERK(i) = trapz(rasmut_pERK_cell{i});
    rasmut_tint_MEK(i)  = trapz(rasmut_MEK_cell{i});
end
figure
plot(max_MEK,max_pERK);
xlabel('Max. MEK');
ylabel('Max. pERK (Molec/cell)');
% High RAF
high_raf_yinit     = yinit;
high_raf_yinit(15) = high_raf_yinit(15) * 10;

highRAF_pERK_cell  = cell(length(k_MEKi),1);
highRAF_MEK_cell   = cell(length(k_MEKi),1);
rasmut_highRAF_pERK_cell  = cell(length(k_MEKi),1);
rasmut_highRAF_MEK_cell   = cell(length(k_MEKi),1);

highRAF_max_pERK   = zeros(length(k_MEKi),1);
highRAF_max_MEK    = zeros(length(k_MEKi),1);

rasmut_highRAF_max_pERK   = zeros(length(k_MEKi),1);
rasmut_highRAF_max_MEK    = zeros(length(k_MEKi),1);

tint_highRAF_pERK  = zeros(length(k_MEKi),1);
tint_highRAF_MEK   = zeros(length(k_MEKi),1);
rasmut_tint_highRAF_pERK  = zeros(length(k_MEKi),1);
rasmut_tint_highRAF_MEK   = zeros(length(k_MEKi),1);
for i=1:length(k_MEKi)
    [~,~,~,~,~,highRAFallValues_wt]             = fullEGFR9_onemodel_MEKi(timeSpan, high_raf_yinit, params, 'min_sor', 'no', k_MEKi(i));
    [~,~,~,~,~,highRAFallValues_rasmut]         = fullEGFR9_onemodel_MEKi(timeSpan, high_raf_yinit, params_mutant, 'min_sor', 'no', k_MEKi(i));
    
    highRAF_pERK_cell{i} = highRAFallValues_wt(:,45); 
    highRAF_MEK_cell{i} = highRAFallValues_wt(:,28); 
    
    rasmut_highRAF_pERK_cell{i} = highRAFallValues_rasmut(:,45); 
    rasmut_highRAF_MEK_cell{i} = highRAFallValues_rasmut(:,28); 
    
    highRAF_max_pERK(i) = max(highRAF_pERK_cell{i});
    highRAF_max_MEK(i)  = max(highRAF_MEK_cell{i});
    
    rasmut_highRAF_max_pERK(i) = max(rasmut_highRAF_pERK_cell{i});
    rasmut_highRAF_max_MEK(i)  = max(rasmut_highRAF_MEK_cell{i});
    
    tint_highRAF_pERK(i) = trapz(highRAF_pERK_cell{i});
    tint_highRAF_MEK(i)  = trapz(highRAF_MEK_cell{i});
    
    rasmut_tint_highRAF_pERK(i) = trapz(rasmut_highRAF_pERK_cell{i});
    rasmut_tint_highRAF_MEK(i)  = trapz(rasmut_highRAF_MEK_cell{i});
end
figure
plot(max_MEK,max_pERK,'LineWidth',2);
hold on
plot(highRAF_max_MEK,highRAF_max_pERK,'LineWidth',2);
xlabel('Max. [MEK] (Molec/cell)');
ylabel('Max. [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');



figure
plot(k_MEKi,max_pERK./max(max_pERK),'LineWidth',2);
hold on
plot(k_MEKi,highRAF_max_pERK./max(highRAF_max_pERK),'LineWidth',2);
xlabel('kMEKi');
ylabel('Normalized Max. [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');

figure
plot(k_MEKi,max_pERK,'LineWidth',2);
hold on
plot(k_MEKi,highRAF_max_pERK,'LineWidth',2);
xlabel('% MEK inhibition');
ylabel('Max. [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');

figure
plot(k_MEKi,tint_pERK./max(tint_pERK),'LineWidth',2);
hold on
plot(k_MEKi,tint_highRAF_pERK./max(tint_highRAF_pERK),'LineWidth',2);
xlabel('k_MEKi');
ylabel('Normalized time-integrated [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');
set(gca, 'xdir', 'reverse');


figure
plot(k_MEKi,rasmut_max_pERK./max(rasmut_max_pERK),'LineWidth',2);
hold on
plot(k_MEKi,rasmut_highRAF_max_pERK./max(rasmut_highRAF_max_pERK),'LineWidth',2);
xlabel('% MEK inhibition');
ylabel('Normalized Max. [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');
title('RAS mutant')

figure
plot(k_MEKi,rasmut_tint_pERK./max(rasmut_tint_pERK),'LineWidth',2);
hold on
plot(k_MEKi,rasmut_tint_highRAF_pERK./max(rasmut_tint_highRAF_pERK),'LineWidth',2);
xlabel('k_MEKi');
ylabel('Normalized time-integrated [pERK] (Molec/cell)');
legend('Low RAF', 'High RAF');
title('RAS mutant')
set(gca, 'xdir', 'reverse');

figure
plot(k_MEKi,max_MEK,'LineWidth',2);
hold on
plot(k_MEKi,highRAF_max_MEK,'LineWidth',2);
xlabel('kMEKi');
ylabel('Max. MEK (Molec/cell)');
legend('Low RAF', 'High RAF');

figure
for i=1:length(k_MEKi)
    txt = ['kMEKi = ',num2str(k_MEKi(i))];
    plot(T, MEK_cell{i},'LineWidth',2,'DisplayName',txt);    
    ylabel('MEK');
    hold on
end
hold off
legend show



figure
for i=1:length(k_MEKi)
    txt = ['kMEKi = ',num2str(k_MEKi(i))];
    plot(T, MEK_cell{i},'LineWidth',2,'DisplayName',txt);  
    plot(T, pERK_cell{i},'LineWidth',2,'DisplayName',txt);    
    ylabel('MEK or pERK');
    xlabel('Time (min)')
    hold on
end
hold off

legend show
