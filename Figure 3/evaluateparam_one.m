function [all_outputs] = evaluateparam_one(wanted_param, params, timeSpan, yinit, number, type)
%function paramOut = evaluateparam_one(wanted_param, params, timeSpan, yinit, number, type)

%rng(1,'twister');
rng(1);

my_other_params_minsor = zeros(number,length(wanted_param));
my_other_params_plussor = zeros(number,length(wanted_param));
my_total_change_params = zeros(number,length(params));
my_other_params_yinit = zeros(1,length(wanted_param)-1);
output = cell(9,1);

my_y_max_plus_RasGTP = zeros(number,1);
my_y_max_plus_Raf1_pm = zeros(number,1);
my_y_max_min_RasGTP = zeros(number,1);
my_y_max_min_Raf1_pm = zeros(number,1);
my_y_max_min_pMEK = zeros(number,1);
my_y_max_min_pERK = zeros(number,1);
my_y_max_min_ERK = zeros(number,1);
my_y_max_min_pEGFR = zeros(number,1);
my_y_max_plussor_pEGFR = zeros(number,1);

my_y_ave_plus_RasGTP = zeros(number,1);
my_y_ave_plus_Raf1_pm = zeros(number,1);
my_y_ave_min_RasGTP = zeros(number,1);
my_y_ave_min_Raf1_pm = zeros(number,1);
my_y_ave_min_pMEK = zeros(number,1);
my_y_ave_min_pERK = zeros(number,1);

my_y_steady_plus_RasGTP = zeros(number,1);
my_y_steady_plus_Raf1_pm = zeros(number,1);
my_y_steady_min_RasGTP = zeros(number,1);
my_y_steady_min_Raf1_pm = zeros(number,1);
my_y_steady_min_pMEK = zeros(number,1);
my_y_steady_min_pERK = zeros(number,1);
my_y_steady_min_ERK = zeros(number,1);
%my_total_steady_plus_RasGTP = zeros(number,1);
t_int_plus_sor_RasGTP = zeros(number,1);
t_int_plus_sor_Raf1_pm = zeros(number,1);
t_int_min_sor_RasGTP = zeros(number,1);
t_int_min_sor_Raf1_pm = zeros(number,1);
t_int_min_sor_pERK = zeros(number,1);
t_int_min_sor_pMEK = zeros(number,1);
t_int_min_sor_ERK = zeros(number,1);
%{
my_y_plus_sor_RasGTP = zeros(number,length(timeSpan));
my_y_plus_sor_Raf1_pm = zeros(number,length(timeSpan));
my_y_min_sor_RasGTP = zeros(number,length(timeSpan));
my_y_min_sor_Raf1_pm = zeros(number,length(timeSpan));
my_y_min_sor_pMEK = zeros(number,length(timeSpan));
my_y_min_sor_pERK = zeros(number,length(timeSpan));
%}

my_y_max_total = cell(6,1); %included ERK
my_y_max_total_Y = cell(6,1);
my_y_ave_total = cell(6,1);
my_y_steady_total = cell(6,1);
t_int_total = cell(6,1);

min_sor_total_output = cell(number,1);
plus_sor_total_output = cell(number, 1);

iteration = 0;
%figure
for k = 1:number
    clc
    iteration = iteration + 1
    if strcmp(type, 'lograndom') == 1
        super = -1 + 2 * rand(length(wanted_param),1);
        super = super';
        changed_params = params';
        y_init_param = yinit;
        my_changed_params = zeros(1,length(wanted_param));
        for i = 1:length(wanted_param)
            changed_params(:,wanted_param(i)) = params(wanted_param(i)) * 10^(super(i));
            my_changed_params(i) = changed_params(:,wanted_param(i));
            
            my_changed_params_yinit = zeros(1,length(wanted_param)-1);
            %{
            for k=2:length(wanted_param)
                y_init_param([22,37,34,40,15,33,28,6]) = my_changed_params(k);
                my_changed_params_yinit(k-1) = my_changed_params(k);
            end
            %}
        end
    end
    
    if strcmp(type, 'random') == 1
        changed_params = params';
        my_changed_params = zeros(1,length(wanted_param));
        min_param = zeros(length(wanted_param),1);
        max_param = zeros(length(wanted_param),1);
        param_range = zeros(number,length(wanted_param));
        y_init_param = yinit;
        for i = 1:length(wanted_param)
            min_param(i,:) = params(wanted_param(i))*10^(-1);
            max_param(i,:) = params(wanted_param(i))*10^(1);
            param_range(:,i) = (max_param(i,:)-min_param(i,:)).*rand(number,1) + min_param(i,:);
            %param_range(:,i) = (max_param(i,:)-min_param(i,:)).*randn(number,1) + min_param(i,:);
            changed_params(wanted_param(i)) = param_range(k,i);
            my_changed_params(i) = changed_params(wanted_param(i));
        end
        
    end
    
    if strcmp(type, 'LHS') == 1
        changed_params = params';
        total_changed_params = repmat(changed_params,number,1);
        my_changed_params = zeros(1,length(wanted_param));
        lhsmatrix = zeros(number,length(wanted_param));
        y_init_param = yinit;
        for i=1:length(wanted_param)
            lhsmatrix(:,i) = LHS_Call(params(wanted_param(i))*0.1, params(wanted_param(i)), params(wanted_param(i))*10, [], number, 'unif');
            %lhsmatrix(:,i) = LHS_Call([], params(wanted_param(i)), [], params(wanted_param(i)) * 0.1, number, 'norm');
            total_changed_params(k,wanted_param(i)) = lhsmatrix(k,i);
            my_changed_params(:,i) = total_changed_params(k, wanted_param(i));
        end
        
        
        %
        %         my_other_params_minsor = zeros(1,length(wanted_param));
        %         my_other_params_plussor = zeros(1,length(wanted_param));
        %         for k = 1:number
        %             for i = 1:length(wanted_param)
        %                 total_changed_params(k,wanted_param(i)) = lhsmatrix(k,i);
        %                 my_changed_params(:,i) = total_changed_params(k, wanted_param(i));
        %             end
        
        
        
    end
%             %plus sor RasGTP max
%             output{1}(k,:) = allValues_rand2(:,184);
%             new_max_RasGTP = max(output{1}(k,:));
%             my_y_max_plus_RasGTP(k) = new_max_RasGTP;
%             
%             %plus sor RasGTP ave
%             new_ave_plus_RasGTP = 1./(T2(end)-T2(1)).*trapz(output{1}(k,:));
%             my_y_ave_plus_RasGTP(k) = new_ave_plus_RasGTP;
%             %plus_sor RasGTP steadystate
%             my_total_steady_plus_RasGTP = allValues_rand2(end,184);
%             my_y_steady_plus_RasGTP(k) = my_total_steady_plus_RasGTP;
%             
%             %plus sor Raf1_pm max
%             output{2}(k,:) = allValues_rand2(:,157); %Raf1_pm_param
%             new_max_Raf1_pm = max(output{2}(k,:));
%             my_y_max_plus_Raf1_pm(k) = new_max_Raf1_pm;
%             
%             %plus sor Raf1_pm ave
%             new_ave_plus_Raf1_pm = 1./(T2(end)-T2(1)).*trapz(output{2}(k,:));
%             my_y_ave_plus_Raf1_pm(k) = new_ave_plus_Raf1_pm;
%             %plus_sor Raf1_pm steadystate
%             my_y_steady_plus_Raf1_pm(k) = allValues_rand2(end,157);
%             
%             my_y_max_total{1} = my_y_max_plus_RasGTP;
%             my_y_max_total{2} = my_y_max_plus_Raf1_pm;
%             
%             my_y_ave_total{1} = my_y_ave_plus_RasGTP;
%             my_y_ave_total{2} = my_y_ave_plus_Raf1_pm;
%             
%             %min sor RasGTP max
%             output{3}(k,:) = allValues_rand(:,184); %RasGTP param
%             new_max_min_RasGTP = max(output{3}(k,:));
%             my_y_max_min_RasGTP(k) = new_max_min_RasGTP;
%             %min sor RasGTP ave
%             new_ave_min_RasGTP = 1./(T1(end)-T1(1)).*trapz(output{3}(k,:));
%             my_y_ave_min_RasGTP(k) = new_ave_min_RasGTP;
%             %steady state
%             my_y_steady_min_RasGTP(k) = allValues_rand(end,184);
%             
%             %min sor Raf1_pm_param max
%             output{4}(k,:) = allValues_rand(:,157); %Raf1_pm_param
%             new_max_min_Raf1_pm = max(output{4}(k,:));
%             my_y_max_min_Raf1_pm(k) = new_max_min_Raf1_pm;
%             %min sor Raf1_pm_param ave
%             new_ave_min_Raf1_pm = 1./(T1(end)-T1(1)).*trapz(output{4}(k,:));
%             my_y_ave_min_Raf1_pm(k) = new_ave_min_Raf1_pm;
%             my_y_steady_min_Raf1_pm(k) = allValues_rand(end,157);
%             
%             %min sor pMEK max
%             output{5}(k,:) = allValues_rand(:,21); %pMEK
%             new_max_pMEK = max(output{5}(k,:));
%             my_y_max_min_pMEK(k) = new_max_pMEK;
%             %min sor pMEK ave
%             new_ave_min_pMEK = 1./(T1(end)-T1(1)).*trapz(output{5}(k,:));
%             my_y_ave_min_pMEK(k) = new_ave_min_pMEK;
%             my_y_steady_min_pMEK(k) = allValues_rand(end,21);
%             
%             %min sor pERK max
%             output{6}(k,:) = allValues_rand(:,45); %pERK
%             new_max_pERK = max(output{6}(k,:));
%             my_y_max_min_pERK(k) = new_max_pERK;
%             %min sor pERK ave
%             new_ave_min_pERK = 1./(T1(end)-T1(1)).*trapz(output{6}(k,:));
%             my_y_ave_min_pERK(k) = new_ave_min_pERK;
%             my_y_steady_min_pERK(k) = allValues_rand(end,45);
%             
%             
%             my_y_max_total{3} = my_y_max_min_RasGTP;
%             my_y_max_total{4} = my_y_max_min_Raf1_pm;
%             my_y_max_total{5} = my_y_max_min_pMEK;
%             my_y_max_total{6} = my_y_max_min_pERK;
% 
%             
%             my_y_ave_total{3} = my_y_ave_min_RasGTP;
%             my_y_ave_total{4} = my_y_ave_min_Raf1_pm;
%             my_y_ave_total{5} = my_y_ave_min_pMEK;
%             my_y_ave_total{6} = my_y_ave_min_pERK;
%             
%             my_y_steady_total{1} = my_y_steady_plus_RasGTP;
%             my_y_steady_total{2} = my_y_steady_plus_Raf1_pm;
%             my_y_steady_total{3} = my_y_steady_min_RasGTP;
%             my_y_steady_total{4} = my_y_steady_min_Raf1_pm;
%             my_y_steady_total{5} = my_y_steady_min_pMEK;
%             my_y_steady_total{6} = my_y_steady_min_pERK;
%             my_y_steady_total{7} = my_y_steady_min_ERK;
%             
%             
%             my_other_params_minsor(k,:) = my_changed_params;
%             my_other_params_plussor(k,:) = my_changed_params;
%             my_total_change_params(k,:) = changed_params;
%             
%             min_sor_total_output{k} = allValues_rand;
%             plus_sor_total_output{k} = allValues_rand2;
%         end
    
    
    if strcmp(type, 'logLHS') == 1
        changed_params = params';
        super = lhsdesign_modified(number,[-1],[1]);
        super = super';
        my_changed_params = zeros(1,length(wanted_param));
        for i = 1:length(wanted_param)
            changed_params(:,wanted_param(i)) = params(wanted_param(i)) * 10^(super(i));
            my_changed_params(i) = changed_params(:,wanted_param(i));
        end
    end
    
    
    if strcmp(type, 'LHS') ~=1
        [T1, ~, ~, params_minsor, ~, allValues_rand]   = fullEGFR9_onemodel(timeSpan, yinit, changed_params,'min_sor', 'vary_init'); %given yinit used, 'vary_init' only specifies the indices of init conditions
        [T2, ~, ~, params_plussor, ~, allValues_rand2] = fullEGFR9_onemodel(timeSpan, yinit, changed_params,'plus_sor', 'vary_init'); %given yinit used, 'vary_init' only specifies the indices of init conditions
    end
    
    if strcmp(type, 'LHS') == 1
        [T1, ~, ~, params_minsor, ~, allValues_rand] = fullEGFR9_onemodel(timeSpan, yinit, total_changed_params(k,:),'min_sor','vary_init');
        [T2, ~, ~, params_plussor, ~, allValues_rand2] = fullEGFR9_onemodel(timeSpan, yinit, total_changed_params(k,:),'plus_sor','vary_init');
    end
    
    %If testing for sensitivity to expression levels:
    if y_init_param ~= yinit
        [T1,~,~,~,~,allValues_rand]  = fullEGFR9_onemodel(timeSpan, y_init_param, params','min_sor');
        [T2,~,~,~,~,allValues_rand2] = fullEGFR9_onemodel(timeSpan, y_init_param, params','plus_sor');
    end
    
    my_changed_params_minsor = zeros(1,length(wanted_param));
    my_changed_params_plussor = zeros(1,length(wanted_param));
    for i = 1:length(wanted_param)
        my_changed_params_minsor(:,i) = params_minsor(wanted_param(i));
        my_changed_params_plussor(:,i) = params_plussor(wanted_param(i));
    end
 
    %plus sor RasGTP max
    output{1}(k,:) = allValues_rand2(:,184);
    new_max_RasGTP = max(output{1}(k,:));
    my_y_max_plus_RasGTP(k) = new_max_RasGTP;
    t_int_plus_sor_RasGTP(k) = trapz(output{1}(k,:));
    %plus sor RasGTP ave
    new_ave_plus_RasGTP = 1./(T2(end)-T2(1)).*trapz(output{1}(k,:));
    my_y_ave_plus_RasGTP(k) = new_ave_plus_RasGTP;
    %plus_sor RasGTP steadystate
    my_total_steady_plus_RasGTP = allValues_rand2(end,184);
    my_y_steady_plus_RasGTP(k) = my_total_steady_plus_RasGTP;
    
    %plus sor Raf1_pm max
    output{2}(k,:) = allValues_rand2(:,157); %Raf1_pm_param
    new_max_Raf1_pm = max(output{2}(k,:));
    my_y_max_plus_Raf1_pm(k) = new_max_Raf1_pm;
    t_int_plus_sor_Raf1_pm(k) = trapz(output{2}(k,:));
    %plus sor Raf1_pm ave
    new_ave_plus_Raf1_pm = 1./(T2(end)-T2(1)).*trapz(output{2}(k,:));
    my_y_ave_plus_Raf1_pm(k) = new_ave_plus_Raf1_pm;
    %plus_sor Raf1_pm steadystate
    my_y_steady_plus_Raf1_pm(k) = allValues_rand2(end,157);
    
    %plus sor pEGFR max
    output{3}(k,:) = allValues_rand2(:,41); %Raf1_pm_param
    new_max_plussor_pEGFR = max(output{3}(k,:));
    my_y_max_plussor_pEGFR(k) = new_max_plussor_pEGFR;

    my_y_max_total{1} = my_y_max_plus_RasGTP;
    my_y_max_total{2} = my_y_max_plus_Raf1_pm;
    %my_y_max_total{3} = my_y_max_plussor_pEGFR;
    
    my_y_ave_total{1} = my_y_ave_plus_RasGTP;
    my_y_ave_total{2} = my_y_ave_plus_Raf1_pm;
   
    %min sor RasGTP max
    output{4}(k,:) = allValues_rand(:,184); %RasGTP param
    new_max_min_RasGTP = max(output{4}(k,:));
    my_y_max_min_RasGTP(k) = new_max_min_RasGTP;
    t_int_min_sor_RasGTP(k) = trapz(output{4}(k,:));
    %min sor RasGTP ave
    new_ave_min_RasGTP = 1./(T1(end)-T1(1)).*trapz(output{4}(k,:));
    my_y_ave_min_RasGTP(k) = new_ave_min_RasGTP;
    %steady state
    my_y_steady_min_RasGTP(k) = allValues_rand(end,184);
    
    %min sor Raf1_pm_param max
    output{5}(k,:) = allValues_rand(:,157); %Raf1_pm_param
    new_max_min_Raf1_pm = max(output{5}(k,:));
    my_y_max_min_Raf1_pm(k) = new_max_min_Raf1_pm;
    %min sor Raf1_pm_param ave
    new_ave_min_Raf1_pm = 1./(T1(end)-T1(1)).*trapz(output{5}(k,:));
    my_y_ave_min_Raf1_pm(k) = new_ave_min_Raf1_pm;
    my_y_steady_min_Raf1_pm(k) = allValues_rand(end,157);
    t_int_min_sor_Raf1_pm(k) = trapz(output{5}(k,:));

    %min sor pMEK max
    output{6}(k,:) = allValues_rand(:,21); %pMEK
    new_max_pMEK = max(output{6}(k,:));
    my_y_max_min_pMEK(k) = new_max_pMEK;
    
    %min sor pMEK ave
    new_ave_min_pMEK = 1./(T1(end)-T1(1)).*trapz(output{5}(k,:));
    my_y_ave_min_pMEK(k) = new_ave_min_pMEK;
    my_y_steady_min_pMEK(k) = allValues_rand(end,21);
    t_int_min_sor_pMEK(k) = trapz(output{6}(k,:));

    %min sor pERK max
    output{7}(k,:) = allValues_rand(:,45); %pERK
    new_max_pERK = max(output{7}(k,:));
    my_y_max_min_pERK(k) = new_max_pERK;
    %min sor pERK ave
    new_ave_min_pERK = 1./(T1(end)-T1(1)).*trapz(output{7}(k,:));
    my_y_ave_min_pERK(k) = new_ave_min_pERK;
    my_y_steady_min_pERK(k) = allValues_rand(end,45);
    t_int_min_sor_pERK(k) = trapz(output{7}(k,:));

    
    %min sor ERK max
    output{8}(k,:) = allValues_rand(:,6); %ERK
    new_max_ERK = max(output{8}(k,:));
    my_y_max_min_ERK(k) = new_max_ERK;
    my_y_steady_min_ERK(k) = allValues_rand(end,6);
    t_int_min_sor_ERK(k) = trapz(output{8}(k,:));
    
    %min sor pEGFR max
    output{9}(k,:) = allValues_rand(:,41); %pEGFR
    new_max_pEGFR = max(output{9}(k,:));
    my_y_max_min_pEGFR(k) = new_max_pEGFR;
    my_y_steady_min_pEGFR(k) = allValues_rand(end,41);
    t_int_min_sor_pEGFR(k) = trapz(output{9}(k,:));

    my_y_max_total{3} = my_y_max_min_RasGTP;
    my_y_max_total{4} = my_y_max_min_Raf1_pm;
    my_y_max_total{5} = my_y_max_min_pMEK;
    my_y_max_total{6} = my_y_max_min_pERK;
    %my_y_max_total{8} = my_y_max_min_pEGFR;
    
    my_y_ave_total{3} = my_y_ave_min_RasGTP;
    my_y_ave_total{4} = my_y_ave_min_Raf1_pm;
    my_y_ave_total{5} = my_y_ave_min_pMEK;
    my_y_ave_total{6} = my_y_ave_min_pERK;
    
    my_y_steady_total{1} = my_y_steady_plus_RasGTP;
    my_y_steady_total{2} = my_y_steady_plus_Raf1_pm;
    my_y_steady_total{3} = my_y_steady_min_RasGTP;
    my_y_steady_total{4} = my_y_steady_min_Raf1_pm;
    my_y_steady_total{5} = my_y_steady_min_pMEK;
    my_y_steady_total{6} = my_y_steady_min_pERK;
    %my_y_steady_total{7} = my_y_steady_min_ERK;

    
    my_other_params_minsor(k,:) = my_changed_params_minsor;
    my_other_params_plussor(k,:) = my_changed_params_plussor;
    my_total_change_params(k,:) = changed_params;
    %my_other_params_yinit(count,:) = my_changed_params_yinit;
    
    min_sor_total_output{k} = allValues_rand;
    plus_sor_total_output{k} = allValues_rand2;
    
    t_int_total{1} = t_int_plus_sor_RasGTP;
    t_int_total{2} = t_int_plus_sor_Raf1_pm;
    t_int_total{3} = t_int_min_sor_RasGTP;
    t_int_total{4} = t_int_min_sor_Raf1_pm;
    t_int_total{5} = t_int_min_sor_pMEK;
    t_int_total{6} = t_int_min_sor_pERK;
    %t_int_total{7} = t_int_min_sor_ERK;
    
    my_y_max_total_Y{1} = my_y_max_plus_RasGTP;
    my_y_max_total_Y{2} = my_y_max_plus_Raf1_pm;
    my_y_max_total_Y{3} = my_y_max_min_RasGTP;
    my_y_max_total_Y{4} = my_y_max_min_Raf1_pm;
    my_y_max_total_Y{5} = my_y_max_min_pMEK;
    my_y_max_total_Y{6} = my_y_max_min_pERK;
    
end

all_outputs.my_other_params_minsor  = my_other_params_minsor;
all_outputs.my_other_params_plussor = my_other_params_plussor;
all_outputs.my_other_params_yinit   = my_other_params_yinit;
all_outputs.my_total_change_params  = my_total_change_params;
all_outputs.my_y_max_total          = my_y_max_total;
all_outputs.my_y_max_total_Y        = my_y_max_total_Y;
all_outputs.my_y_ave_total          = my_y_ave_total;
all_outputs.my_y_steady_total       = my_y_steady_total;
all_outputs.t_int_total             = t_int_total;
all_outputs.output                  = output;
all_outputs.min_sor_total_output    = min_sor_total_output;
all_outputs.plus_sor_total_output   = plus_sor_total_output;

end



