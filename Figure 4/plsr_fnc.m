function [plsr_z_score_x, plsr_z_score_y, Xloadings, vipScores, Q2Y, R2Y, BETA, PCTVAR, PC1loadings] = plsr_fnc(number, type, raw_x, raw_y_cell, ncomp, wanted_param)


%Pre-process independent and dependent variables used in PLSR 
%1) If 'log-random' or 'log-LHS' sampling is used, input values
%are log-transformed then z-scored.
%2) If 'random' or 'LHS' sampling is used, input values are z-scored.
%Inputs: 
    % number: number of parameter sets
    % type: sampling method used to generate parameter sets
    % raw_x: parameter sets
    % raw_y_cell: ODE solutions
    % ncomp: number of PLS components
    % wanted_param: parameters of interest (rate constants)
    
%Outputs: 
    % plsr_z_score_x: z-scored parameter values 
    % plsr_z_score_y: z-scored ODE solutions 
    % Xloadings: Parameter loadings in each PLS component from PLSR SA 
    % vipScores: VIP scores from PLSR SA
    % Q2Y: Predictive power of PLSR SA model using leave-one-out cross
    % validation 
    % R2Y: Explanatory power of PLSR SA model
    % BETA: Coefficient estimates of each parameter in PLSR SA 
    % PCTVAR: Percent variance explained by PLSR SA model
    % PC1loadings: Parameter loadings in the first PLS component, normalized by the maximum loading  



%Convert cell to double if necessary 
raw_varied_params = raw_x;
size(raw_varied_params,2)
if isa(raw_y_cell,'cell')
    raw_y = [raw_y_cell{:}];
else
    raw_y = raw_y_cell;
end


if strcmp(type, 'lograndom') == 1
    param_size = size(raw_varied_params);
    observations = param_size(1);
    my_log_params = zeros(param_size);
    for o = 1:observations
        for p = 1:size(raw_varied_params,2)
            if raw_varied_params(o,p) == 0
                my_log_params(o,p) = raw_varied_params(o,p);
            else
                my_log_params(o,p) = log10(raw_varied_params(o,p));
            end
        end
    end
    
    plsr_z_score_x = zeros(number, size(raw_varied_params,2)); 
    for p = 1:size(raw_varied_params,2)
        x = my_log_params(:,p);
        Zsc = @(x) (x - mean(x))./std(x);  
        Zx = Zsc(x)  ;
        plsr_z_score_x(:,p) = Zx;
    end
    
    plsr_z_score_y = zeros(number, size(raw_y,2)); 
    for p = 1:size(raw_y,2)
        y = raw_y(:,p);
        Zsc = @(y) (y - mean(y))./std(y);   
        Zy = Zsc(y)  ;
        plsr_z_score_y(:,p) = Zy;
    end
end

if strcmp(type, 'logLHS') == 1
    param_size = size(raw_varied_params);
    observations = param_size(1);
    my_log_params = zeros(param_size);
    for o = 1:observations
        for p = 1:size(raw_varied_params,2)
            if raw_varied_params(o,p) == 0
                my_log_params(o,p) = raw_varied_params(o,p);
            else
                my_log_params(o,p) = log10(raw_varied_params(o,p));
            end
        end
    end
    plsr_z_score_x = zeros(number, size(raw_varied_params,2)); 
    for p = 1:size(raw_varied_params,2)
        x = my_log_params(:,p);
        Zsc = @(x) (x - mean(x))./std(x);   
        Zx = Zsc(x)  ;
        plsr_z_score_x(:,p) = Zx;
    end
    plsr_z_score_y = zeros(number, size(raw_y,2)); 
    for p = 1:size(raw_y,2)
        y = raw_y(:,p);
        Zsc = @(y) (y - mean(y))./std(y);   
        Zy = Zsc(y)  ;
        plsr_z_score_y(:,p) = Zy;
    end
end

if strcmp(type, 'random') == 1
    plsr_z_score_x = zeros(number,size(raw_varied_params,2)); 
    for p = 1:size(raw_varied_params,2)
        x = raw_varied_params(:,p);
        Zsc = @(x) (x - mean(x))./std(x);   
        Zx = Zsc(x)  ;
        plsr_z_score_x(:,p) = Zx;
    end
    plsr_z_score_y = zeros(number, size(raw_y,2)); 
    for k = 1:size(raw_y,2)
        y = raw_y(:,k);
        Zsc = @(y) (y - mean(y))./std(y);   
        Zy = Zsc(y);
        plsr_z_score_y(:,k) = Zy;
    end
end
if strcmp(type, 'LHS') == 1
    plsr_z_score_x = zeros(number, size(raw_varied_params,2)); 
    for p = 1:size(raw_varied_params,2)
        x = raw_varied_params(:,p);
        Zsc = @(x) (x - mean(x))./std(x);   
        Zx = Zsc(x)  ;
        plsr_z_score_x(:,p) = Zx;
    end
    plsr_z_score_y = zeros(number, size(raw_y,2)); 
    for p = 1:size(raw_y,2)
        y = raw_y(:,p);
        Zsc = @(y) (y - mean(y))./std(y);  
        Zy = Zsc(y)  ;
        plsr_z_score_y(:,p) = Zy;
    end
end


C1 = cvpartition(number, 'LeaveOut');
[Xloadings,YL,XS,YS,BETA,PCTVAR,PLSmsep,stats] = plsregress(plsr_z_score_x, plsr_z_score_y, ncomp,'CV',C1);
Q2Y = 1- PLSmsep(2,2:end)/sum(sum((plsr_z_score_y-mean(plsr_z_score_y)).^2)./size(plsr_z_score_y,1));
R2Y = cumsum(PCTVAR(2,1:end));
R2X = cumsum(PCTVAR(1,1:end));
yfitPLS = [ones(number,1) plsr_z_score_x]*BETA;
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
vipScores = sqrt(size(Xloadings,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));


max_loadings = max(abs(Xloadings(:,1)));
sensitivities = zeros(size(Xloadings(:,1)));
for p = 1:length(wanted_param)
    sensitivities(p) = abs(Xloadings(p,1)) / max_loadings * 100;
end

PC1loadings = sensitivities;

end

