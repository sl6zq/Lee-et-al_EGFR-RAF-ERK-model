% Datasets from Pinilla-Macua et al., 2016 (https://doi.org/10.1073/pnas.1520301113) 
% Surve et al., 2019 (https://www.molbiolcell.org/doi/full/10.1091/mbc.E18-08-0512) used for model fitting
%% Define the number of data points to generate (for CI calculation in Fig.5)
datasample = 20; 
rng(0); %random seed for synthetic data generation (for CI generation in Fig. 5)

%% RAS-GTP (EGF 10ng/mL), Pinilla-Macua et al., 2016
min_sor_rastimedata                 =  [2.5 5.0 10.0 15.0 30.0 60.0]';
min_sor_rasdatapcnts                = [99.46 94.34 61.31 13.59 5.53 14.51]';
min_sor_rasdatanums                 = min_sor_rasdatapcnts.*270;
min_sor_rasdataerror                = [1917 1971 4185 1944 1053 2349]';
min_sor_rasdataerrorpercent         = min_sor_rasdataerror./270;
min_sor_rasdata_rangenums           = zeros(datasample,length(min_sor_rastimedata));
for a = 1:length(min_sor_rastimedata)
    min_sor_rasdata_rangenums(:,a)  = normrnd(min_sor_rasdatanums(a),min_sor_rasdataerror(a),[datasample 1]);
end
%% Membrane RAF1 (EGF 10ng/mL), Surve et al., 2019
min_sor_raftimedata                 = [0.5 0.75 1.5 2.25 2.5 3 3.5 4 5 5.5 6 7 7.5 9.5 10 10.5 14 16 17 19 21 23 24 26 29]';
min_sor_rafdatanums                 = [-24 668 341 1026 1219 1261 869 1730 999 526 1771 217 77 3 337 99 80 -88.2 38 -33 71 202 213 161 106]';
min_sor_rafdataerror                = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200]';  
min_sor_rafdata_rangenums           = zeros(datasample,length(min_sor_raftimedata));
for b = 1:length(min_sor_raftimedata)
    min_sor_rafdata_rangenums(:,b)  = normrnd(min_sor_rafdatanums(b),min_sor_rafdataerror(b),[datasample 1]);
end
%% pERK (EGF 10ng/mL), Surve et al., 2019
min_sor_perktimedata                = [5 10 20 40 60]';
min_sor_perkdatapcnts               = [95.4 88.2 77.2 56.2 51.5]';
min_sor_perkdatanums                = min_sor_perkdatapcnts.*2100;
min_sor_perkdataerrorpercent        = [0.0682 0.113 0.141 0.265 0.235]';
min_sor_perkdataerror               = [0.0682 0.113 0.141 0.265 0.235]' .* 210000;
min_sor_perkdata_rangenums          = zeros(datasample,length(min_sor_perktimedata));
for c = 1:length(min_sor_perktimedata)
    min_sor_perkdata_rangenums(:,c) = normrnd(min_sor_perkdatanums(c),min_sor_perkdataerror(c),[datasample 1]);
end
%% pMEK (EGF 10ng/mL), Surve et al., 2019
min_sor_pmektimedata                = [5 10 20 40 60]';
min_sor_pmekdatapcnts               = [95.6 74.6 46.8 54.9 74.5]'; 
min_sor_pmekdatanums                = min_sor_pmekdatapcnts.*1750;
min_sor_pmekdataerror               = [0.117 0.194 0.157 0.0685 0.0596]' .* 175000;
min_sor_pmekdataerrorpercent        = [0.117 0.194 0.157 0.0685 0.0596]';
min_sor_pmekdata_rangenums          = zeros(datasample,length(min_sor_pmektimedata));
for c = 1:length(min_sor_pmektimedata)
    min_sor_pmekdata_rangenums(:,c) = normrnd(min_sor_pmekdatanums(c),min_sor_pmekdataerror(c),[datasample 1]);
end
%% RAS-GTP (EGF 10ng/mL, 10 μM sorafenib), assumed - see discussion in Surve et al., 2019
plus_sor_rastimedata                = [2.5 5.0 10.0 15.0 30.0 60.0]';
plus_sor_rasdatapcnts               = [99.46 94.34 61.31 13.59 5.53 14.51]';
plus_sor_rasdatanums                = plus_sor_rasdatapcnts.*270;
plus_sor_rasdataerror               = [1917 1971 4185 1944 1053 2349]';
plus_sor_rasdataerrorpercent        = plus_sor_rasdataerror./270;
plus_sor_rasdata_rangenums          = zeros(datasample,length(plus_sor_rastimedata));
for d = 1:length(plus_sor_rastimedata)
    plus_sor_rasdata_rangenums(:,d) = normrnd(plus_sor_rasdatanums(d),plus_sor_rasdataerror(d),[datasample 1]);
end
%% Membrane RAF1 (EGF 10ng/mL, 10 μM sorafenib), Surve et al., 2019
plus_sor_raftimedata                = [0.5 2 4 6 7 9 9.5 12 15 15.5 17 19 20 21.5 24 28 28.5]';
plus_sor_rafdatanums                = [356 487 270 2215 1906 3405 1548 2466 2052 2939 3509 1732 3594 2219 2061 2258 1846]';
plus_sor_rafdataerror               = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200]';  
plus_sor_rafdataerrorpercent        = plus_sor_rafdataerror./11.7;
plus_sor_rafdata_rangenums          = zeros(datasample,length(plus_sor_raftimedata));
for e = 1:length(plus_sor_raftimedata)
    plus_sor_rafdata_rangenums(:,e) = normrnd(plus_sor_rafdatanums(e),plus_sor_rafdataerror(e),[datasample 1]);
end
%% Store mean data points in structure
dat_struc.timepoints{1} = min_sor_rastimedata;
dat_struc.timepoints{2} = min_sor_raftimedata;
dat_struc.timepoints{3} = min_sor_perktimedata;
dat_struc.timepoints{4} = min_sor_pmektimedata;
dat_struc.timepoints{5} = plus_sor_rastimedata;
dat_struc.timepoints{6} = plus_sor_raftimedata;

dat_struc.datapoints{1} = min_sor_rasdatanums;
dat_struc.datapoints{2} = min_sor_rafdatanums;
dat_struc.datapoints{3} = min_sor_perkdatanums;
dat_struc.datapoints{4} = min_sor_pmekdatanums;
dat_struc.datapoints{5} = plus_sor_rasdatanums;
dat_struc.datapoints{6} = plus_sor_rafdatanums;


dat_struc.dataerror{1} = min_sor_rasdataerror;
dat_struc.dataerror{2} = min_sor_rafdataerror;
dat_struc.dataerror{3} = min_sor_perkdataerror;
dat_struc.dataerror{4} = min_sor_pmekdataerror;
dat_struc.dataerror{5} = plus_sor_rasdataerror;
dat_struc.dataerror{6} = plus_sor_rafdataerror;


