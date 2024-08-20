function [T,Y,yinit,param, allNames, allValues] = fullEGFR9_onemodel(argTimeSpan,argYinit,argParam, modeltype, yinit_input)

% 2.17.21

% [T,Y,yinit,param] = fullEGFR9_onemodel(argTimeSpan,argYinit,argParam)

% input:
%     argTimeSpan is a vector of start and stop times (e.g. timeSpan = [0 10.0])
%     argYinit is a vector of initial conditions for the state variables (optional)
%     argParam is a vector of values for the parameters (optional)
%
% output:
%     T is the vector of times
%     Y is the vector of state variables
%     yinit is the initial conditions that were used
%     param is the parameter vector that was used
%     allNames is the output solution variable names
%     allValues is the output solution variable values corresponding to the names
%
%     example of running this file: [T,Y,yinit,param,allNames,allValues] = myMatlabFunc; <-(your main function name)


% output variable lengh and names
numVars = 191;
allNames = {'Ras_iBRaf';'Ras_Raf1_iBRaf_tetramer';'Ras_BRaf';'Ras_nfpRaf1';'iRaf1';'ERK';'nfpSOS';'Ras_pRaf1';'BRaf_iBRaf_dimer';'Ras_pRaf1_iBRaf_tetramer';'Ras_pRaf1_tetramer';'Ras_nfpiRaf1';'Ras_GTP';'nfpiRaf1';'Raf1';'Ras_BRaf_pRaf1_tetramer';'Ras_pRaf1_Raf1_tetramer';'nfpRaf1';'Ras_nfpBRaf';'iBRaf';'pMEK';'mE';'mEL';'pRaf1';'EG2';'Ras_iRaf1';'Ras_nfpiBRaf';'MEK';'Ras_Braf_iRaf1_tetramer';'Ras_Raf1';'mELmEL';'nfpiBRaf';'BRaf';'SOS';'nfpBRaf';'GRB2_SOS';'GRB2';'Ras_iRaf1_tetramer';'iBRaf_dimer';'Ras_GDP';'E';'Ras_BRaf_Raf1_tetramer';'BRaf_dimer';'EG2SOS';'pERK';'Ras_iRaf1_iBRaf1_tetramer';'Kr_E_G2_bind';'Ras_GDP_init_molecules_um_2';'Kf_iRaf_iBRaf_dim';'Raf1_cyt_param';'Kr_E_G2SOS_bind';'Kf_E_G2SOS_bind';'J_E_G2SOS_bind';'Kr_BRAF_dim';'Kf_RAS_GTP_bind';'J_RAS_GTP_bind';'Kf_Ras_iRaf_nfp';'J_Ras_iRaf_nfp';'Kf_Ras_Raf1_nfp';'Kf_Ras_iBRaf_nfp';'Kf_pRAF1_dephos';'J_pRAF1_dephos';'Kf_iRaf1_iRaf1_bind';'Kf_BRAF_RAF1_phos';'Kf_BRAF_dim';'J_BRAF_dim';'sorafenib_init_molecules_um_3';'sorafenib';'Kf_E_G2_bind';'J_E_G2_bind';'Kf_BRAF_nfp';'Kf_neg_fdbk_dephos_pSOS';'J_neg_fdbk_dephos_pSOS';'Kr_Ras_Raf1_bind';'Kf_neg_fdbk_dephos_iBRAF';'Kf_ERK_dephos';'Kf_iBRaf_dim';'pRaf1_total';'Kf_EGF_bind';'Kf_EG2_SOS_bind';'Kf_RAS_GTPhydro';'J_RAS_GTPhydro';'Kf_iBRaf_Raf1_bind';'Kf_neg_fdbk_unbind_RASt_nfpBRAF';'Kr_mELmEL_phos_dephos';'Kf_RAS_BRAF_bind';'Kr_iBRaf_Raf1_bind';'Kf_pRaf1_Raf1_bind';'Kf_neg_fdbk_phos_of_SOS';'Kr_RAS_pRAF1_unbind';'Kf_neg_fdbk_unbind_piRaf1_RASt';'Kf_Ras_Raf1_bind';'J_Ras_Raf1_bind';'Kf_neg_fdbk_unbind_piBRaf_RASt';'Kf_Raf1_iBRaf_phos';'J_Raf1_iBRaf_phos';'Kr_BRaf_iBRaf_dim';'Kf_BRaf_iBRaf_dim';'Kf_sorafenib_Raf1_bind_cyt';'J_ERK_dephos';'Kf_Ras_pRaf1_dp';'Kf_GRB2_SOS_bind_cyt';'Kr_GRB2_SOS_bind_cyt';'J_GRB2_SOS_bind_cyt';'Kf_mEL_dim';'Kf_MEK_dephos';'J_neg_fdbk_unbind_RASt_nfpBRAF';'Kr_iRaf_iBRaf_dim';'Kr_RAS_iBRAF_bind';'Kf_RAS_iBRAF_bind';'J_RAS_iBRAF_bind';'Kf_neg_fdbk_dephos_pBRAF';'J_neg_fdbk_dephos_pBRAF';'Kf_Ras_iRaf1_bind';'Kr_Ras_iRaf1_bind';'J_Ras_iRaf1_bind';'SOS_init_molecules_um_3';'mE_init_molecules_um_2';'KFlux_pm_cytoplasm';'ERK_init_molecules_um_3';'Kr_EG2_SOS_bind';'J_EG2_SOS_bind';'J_MEK_dephos';'Kf_ERK_phos';'J_ERK_phos';'Kf_BRAF_RAF1_bind';'Kr_BRAF_RAF1_bind';'J_BRAF_RAF1_bind';'Kf_pRaf1_pRaf1_unbind';'Kr_pRaf1_pRaf1_unbind';'J_pRaf1_pRaf1_unbind';'active_Raf';'Kf_MEK_phos';'J_MEK_phos';'Kf_BRaf_pRaf1_unbind';'Kr_BRaf_pRaf1_unbind';'J_BRaf_pRaf1_unbind';'Kf_BRaf_iRaf1_bind';'Kf_pRaf1_iBRaf_unbind';'Kr_pRaf1_iBRaf_unbind';'J_pRaf1_iBRaf_unbind';'Kf_sorafenib_BRaf_bind_cyt';'J_iBRaf_Raf1_bind';'Kr_mEL_dim';'J_mEL_dim';'Kr_pRaf1_Raf1_bind';'Kf_pRAF1_RAF1_phos';'Kf_RAS_pRAF1_unbind';'J_RAS_pRAF1_unbind';'Kr_BRaf_iRaf1_bind';'Kf_neg_fdbck_dephos_piRaf';'J_Ras_Raf1_nfp';'Kf_neg_fdbk_dephos_pRaf1';'J_neg_fdbk_dephos_pRaf1';'Kr_sorafenib_Raf1_bind_cyt';'J_sorafenib_Raf1_bind_cyt';'Raf1_pm_param';'boundRaf_frac_parameter';'Kr_EGF_bind';'J_BRAF_nfp';'J_Ras_iBRaf_nfp';'Kr_RAS_BRAF_bind';'J_EGF_bind';'Kr_iRaf1_iRaf1_bind';'MEK_init_molecules_um_3';'J_iRaf1_iRaf1_bind';'GRB2_init_molecules_um_3';'J_BRaf_iBRaf_dim';'Kf_mELmEL_phos_dephos';'J_mELmEL_phos_dephos';'J_RAS_BRAF_bind';'J_iRaf_iBRaf_dim';'J_neg_fdbk_dephos_iBRAF';'J_Ras_pRaf1_dp';'Kr_iBRaf_dim';'Kf_neg_fdbk_unbind_RASt_nfpRAF1';'Raf1_init_molecules_um_3';'J_BRAF_RAF1_phos';'Kr_sorafenib_BRaf_bind_cyt';'J_sorafenib_BRaf_bind_cyt';'J_pRaf1_Raf1_bind';'J_neg_fdbk_unbind_piRaf1_RASt';'BRaf_init_molecules_um_3';'RasGTP_parameter';'J_BRaf_iRaf1_bind';'J_neg_fdbk_phos_of_SOS';'J_neg_fdbk_unbind_piBRaf_RASt';'J_iBRaf_dim';'J_neg_fdbck_dephos_piRaf';'J_pRAF1_RAF1_phos';'J_neg_fdbk_unbind_RASt_nfpRAF1';};

if nargin >= 3
	if length(argParam) > 0
		%
		% parameter values overridden by function arguments
		%
		param = argParam;
	end
end
if nargin >= 1
	if length(argTimeSpan) > 0
		%
		% TimeSpan overridden by function arguments
		%
		timeSpan = argTimeSpan;
	end
end

if strcmp(yinit_input, 'vary_yinit') == 1
    
    ERK_init_molecules_um_3 = param(19);
    Raf1_init_molecules_um_3 = param(115);
    mE_init_molecules_um_2 = param(141);
    MEK_init_molecules_um_3 = param(31);
    BRaf_init_molecules_um_3 = param(131);
    SOS_init_molecules_um_3 = param(11);
    GRB2_init_molecules_um_3 = param(109);
    Ras_GDP_init_molecules_um_2 = param(138);
    
    yinit = [
        0.0;		% yinit(1) is the initial condition for 'Ras_iBRaf'
        0.0;		% yinit(2) is the initial condition for 'Ras_Raf1_iBRaf_tetramer'
        0.0;		% yinit(3) is the initial condition for 'Ras_BRaf'
        0.0;		% yinit(4) is the initial condition for 'Ras_nfpRaf1'
        0.0;		% yinit(5) is the initial condition for 'iRaf1'
        ERK_init_molecules_um_3;		% yinit(6) is the initial condition for 'ERK'
        0.0;		% yinit(7) is the initial condition for 'nfpSOS'
        0.0;		% yinit(8) is the initial condition for 'Ras_pRaf1'
        0.0;		% yinit(9) is the initial condition for 'BRaf_iBRaf_dimer'
        0.0;		% yinit(10) is the initial condition for 'Ras_pRaf1_iBRaf_tetramer'
        0.0;		% yinit(11) is the initial condition for 'Ras_pRaf1_tetramer'
        0.0;		% yinit(12) is the initial condition for 'Ras_nfpiRaf1'
        0.0;		% yinit(13) is the initial condition for 'Ras_GTP'
        0.0;		% yinit(14) is the initial condition for 'nfpiRaf1'
        Raf1_init_molecules_um_3;		% yinit(15) is the initial condition for 'Raf1'
        0.0;		% yinit(16) is the initial condition for 'Ras_BRaf_pRaf1_tetramer'
        0.0;		% yinit(17) is the initial condition for 'Ras_pRaf1_Raf1_tetramer'
        0.0;		% yinit(18) is the initial condition for 'nfpRaf1'
        0.0;		% yinit(19) is the initial condition for 'Ras_nfpBRaf'
        0.0;		% yinit(20) is the initial condition for 'iBRaf'
        0.0;		% yinit(21) is the initial condition for 'pMEK'
        mE_init_molecules_um_2;		% yinit(22) is the initial condition for 'mE'
        0.0;		% yinit(23) is the initial condition for 'mEL'
        0.0;		% yinit(24) is the initial condition for 'pRaf1'
        0.0;		% yinit(25) is the initial condition for 'EG2'
        0.0;		% yinit(26) is the initial condition for 'Ras_iRaf1'
        0.0;		% yinit(27) is the initial condition for 'Ras_nfpiBRaf'
        MEK_init_molecules_um_3;		% yinit(28) is the initial condition for 'MEK'
        0.0;		% yinit(29) is the initial condition for 'Ras_Braf_iRaf1_tetramer'
        0.0;		% yinit(30) is the initial condition for 'Ras_Raf1'
        0.0;		% yinit(31) is the initial condition for 'mELmEL'
        0.0;		% yinit(32) is the initial condition for 'nfpiBRaf'
        BRaf_init_molecules_um_3;		% yinit(33) is the initial condition for 'BRaf'
        SOS_init_molecules_um_3;		% yinit(34) is the initial condition for 'SOS'
        0.0;		% yinit(35) is the initial condition for 'nfpBRaf'
        0.0;		% yinit(36) is the initial condition for 'GRB2_SOS'
        GRB2_init_molecules_um_3;		% yinit(37) is the initial condition for 'GRB2'
        0.0;		% yinit(38) is the initial condition for 'Ras_iRaf1_tetramer'
        0.0;		% yinit(39) is the initial condition for 'iBRaf_dimer'
        Ras_GDP_init_molecules_um_2;		% yinit(40) is the initial condition for 'Ras_GDP'
        0.0;		% yinit(41) is the initial condition for 'E'
        0.0;		% yinit(42) is the initial condition for 'Ras_BRaf_Raf1_tetramer'
        0.0;		% yinit(43) is the initial condition for 'BRaf_dimer'
        0.0;		% yinit(44) is the initial condition for 'EG2SOS'
        0.0;		% yinit(45) is the initial condition for 'pERK'
        0.0;		% yinit(46) is the initial condition for 'Ras_iRaf1_iBRaf1_tetramer'
        ];
else
    yinit = argYinit;
end


% Default Initial Conditions
%
%{
yinit = [
	0.0;		% yinit(1) is the initial condition for 'Ras_iBRaf'
	0.0;		% yinit(2) is the initial condition for 'Ras_Raf1_iBRaf_tetramer'
	0.0;		% yinit(3) is the initial condition for 'Ras_BRaf'
	0.0;		% yinit(4) is the initial condition for 'Ras_nfpRaf1'
	0.0;		% yinit(5) is the initial condition for 'iRaf1'
	ERK_init_molecules_um_3;		% yinit(6) is the initial condition for 'ERK'
	0.0;		% yinit(7) is the initial condition for 'nfpSOS'
	0.0;		% yinit(8) is the initial condition for 'Ras_pRaf1'
	0.0;		% yinit(9) is the initial condition for 'BRaf_iBRaf_dimer'
	0.0;		% yinit(10) is the initial condition for 'Ras_pRaf1_iBRaf_tetramer'
	0.0;		% yinit(11) is the initial condition for 'Ras_pRaf1_tetramer'
	0.0;		% yinit(12) is the initial condition for 'Ras_nfpiRaf1'
	0.0;		% yinit(13) is the initial condition for 'Ras_GTP'
	0.0;		% yinit(14) is the initial condition for 'nfpiRaf1'
	Raf1_init_molecules_um_3;		% yinit(15) is the initial condition for 'Raf1'
	0.0;		% yinit(16) is the initial condition for 'Ras_BRaf_pRaf1_tetramer'
	0.0;		% yinit(17) is the initial condition for 'Ras_pRaf1_Raf1_tetramer'
	0.0;		% yinit(18) is the initial condition for 'nfpRaf1'
	0.0;		% yinit(19) is the initial condition for 'Ras_nfpBRaf'
	0.0;		% yinit(20) is the initial condition for 'iBRaf'
	0.0;		% yinit(21) is the initial condition for 'pMEK'
	mE_init_molecules_um_2;		% yinit(22) is the initial condition for 'mE'
	0.0;		% yinit(23) is the initial condition for 'mEL'
	0.0;		% yinit(24) is the initial condition for 'pRaf1'
	0.0;		% yinit(25) is the initial condition for 'EG2'
	0.0;		% yinit(26) is the initial condition for 'Ras_iRaf1'
	0.0;		% yinit(27) is the initial condition for 'Ras_nfpiBRaf'
	MEK_init_molecules_um_3;		% yinit(28) is the initial condition for 'MEK'
	0.0;		% yinit(29) is the initial condition for 'Ras_Braf_iRaf1_tetramer'
	0.0;		% yinit(30) is the initial condition for 'Ras_Raf1'
	0.0;		% yinit(31) is the initial condition for 'mELmEL'
	0.0;		% yinit(32) is the initial condition for 'nfpiBRaf'
	BRaf_init_molecules_um_3;		% yinit(33) is the initial condition for 'BRaf'
	SOS_init_molecules_um_3;		% yinit(34) is the initial condition for 'SOS'
	0.0;		% yinit(35) is the initial condition for 'nfpBRaf'
	0.0;		% yinit(36) is the initial condition for 'GRB2_SOS'
	GRB2_init_molecules_um_3;		% yinit(37) is the initial condition for 'GRB2'
	0.0;		% yinit(38) is the initial condition for 'Ras_iRaf1_tetramer'
	0.0;		% yinit(39) is the initial condition for 'iBRaf_dimer'
	Ras_GDP_init_molecules_um_2;		% yinit(40) is the initial condition for 'Ras_GDP'
	0.0;		% yinit(41) is the initial condition for 'E'
	0.0;		% yinit(42) is the initial condition for 'Ras_BRaf_Raf1_tetramer'
	0.0;		% yinit(43) is the initial condition for 'BRaf_dimer'
	0.0;		% yinit(44) is the initial condition for 'EG2SOS'
	0.0;		% yinit(45) is the initial condition for 'pERK'
	0.0;		% yinit(46) is the initial condition for 'Ras_iRaf1_iBRaf1_tetramer'
];
%}
%{
if nargin >= 2
	if length(argYinit) > 0
		%
		% initial conditions overridden by function arguments
		%
		yinit = argYinit;
	end
end
%}
%
% Default Parameters
%   constants are only those "Constants" from the Math Description that are just floating point numbers (no identifiers)
%   note: constants of the form "A_init" are really initial conditions and are treated in "yinit"
%
%{
param = [
	0.0;		% param(1) is 'Ras_iRaf1_tetramer_init_molecules_um_2'
	7.65E+03; %310.0;		% param(2) is 'kEf'
	96485.3321;		% param(3) is 'mlabfix_F_'
	0.0;		% param(4) is 'Ras_pRaf1_Raf1_tetramer_init_molecules_um_2'
	2.24E+01; %6.60624874;		% param(5) is 'kcatE'
	0.0;		% param(6) is 'pRaf1_init_molecules_um_3'
	1.0;		% param(7) is 'netValence_neg_fdbk_phos_of_SOS'
	3.20E-03; %6.0E-6;		% param(8) is 'kpMEK'
	1.0;		% param(9) is 'netValence_E_G2SOS_bind'
	5.06E-05; %0.0025;		% param(10) is 'knfpSOS'
	7700.0;		% param(11) is 'init_SOS_HeLa'
	1.0;		% param(12) is 'netValence_pRaf1_Raf1_bind'
	1.0;		% param(13) is 'netValence_iBRaf_dim'
	1.0;		% param(14) is 'netValence_BRaf_iRaf1_bind'
	1.0;		% param(15) is 'netValence_neg_fdbk_unbind_RASt_nfpRAF1'
	1.65E-04; %0.037;		% param(16) is 'kBr'
	0.0;		% param(17) is 'Kr_ERK_phos'
	0.0;		% param(18) is 'Kr_RAS_GTP_bind'
	600000.0;		% param(19) is 'init_ERK_HeLa'
	0.0;		% param(20) is 'Kr_Ras_iRaf_nfp'
	1.0;		% param(21) is 'netValence_pRAF1_RAF1_phos'
	1.73E-07; %2.5E-7;		% param(22) is 'kBf'
	4.49E+00; %0.6;		% param(23) is 'kdpERK'
	0.0;		% param(24) is 'Kr_pRAF1_dephos'
	1.14E+03; %8.4;		% param(25) is 'kR1r'
	0.0;		% param(26) is 'iRaf1_init_molecules_um_3'
	3.79E-03; %2.5E-5;		% param(27) is 'kR1f'
	1.0;		% param(28) is 'netValence_RAS_pRAF1_unbind'
	0.0;		% param(29) is 'Kr_neg_fdbk_unbind_RASt_nfpBRAF'
	0.0;		% param(30) is 'Kr_MEK_phos'
	500000.0;		% param(31) is 'init_MEK_HeLa'
	1.0;		% param(32) is 'netValence_Ras_iBRaf_nfp'
	0.0;		% param(33) is 'Kr_neg_fdbk_dephos_pRaf1'
	1.0;		% param(34) is 'netValence_RAS_iBRAF_bind'
	0.0;		% param(35) is 'nfpiBRaf_init_molecules_um_3'
	1.0;		% param(36) is 'netValence_E_G2_bind'
	0.0;		% param(37) is 'Kr_RAS_GTPhydro'
	0.0;		% param(38) is 'Kr_Ras_iBRaf_nfp'
	0.0;		% param(39) is 'pMEK_init_molecules_um_3'
	0.0;		% param(40) is 'Ras_BRaf_init_molecules_um_2'
	0.0;		% param(41) is 'Ras_iRaf1_iBRaf1_tetramer_init_molecules_um_2'
	0.0;		% param(42) is 'Ras_pRaf1_tetramer_init_molecules_um_2'
	1000.0;		% param(43) is 'K_millivolts_per_volt'
	1.0E-9;		% param(44) is 'mlabfix_K_GHK_'
	1.68E-16; %6.0E-7;		% param(45) is 'knfpBR1'
	0.0;		% param(46) is 'EG2_init_molecules_um_2'
	0.0;		% param(47) is 'Ras_BRaf_pRaf1_tetramer_init_molecules_um_2'
	1.0;		% param(48) is 'netValence_pRaf1_pRaf1_unbind'
	3.97E+00; %190.0;		% param(49) is 'kpR1'
	7.67E-03; %0.0083;		% param(50) is 'kSOSr'
	0.0;		% param(51) is 'Ras_nfpRaf1_init_molecules_um_2'
	1.0;		% param(52) is 'netValence_iRaf1_iRaf1_bind'
	4.22E-07; %2.2E-6;		% param(53) is 'kSOSf'
	0.0;		% param(54) is 'Ras_iBRaf_init_molecules_um_2'
	1.0;		% param(55) is 'netValence_RAS_BRAF_bind'
	0.0;		% param(56) is 'Kr_ERK_dephos'
	1.38E+01; %5.838226569;		% param(57) is 'kRgneslow'
	9.64853321E-5;		% param(58) is 'mlabfix_F_nmol_'
	1.0;		% param(59) is 'netValence_mEL_dim'
	1.0;		% param(60) is 'netValence_BRAF_RAF1_bind'
	0.0;		% param(61) is 'GRB2_SOS_init_molecules_um_3'
	0.0;		% param(62) is 'Kr_BRAF_RAF1_phos'
	0.0;		% param(63) is 'E_init_molecules_um_2'
	0.0;		% param(64) is 'Ras_GTP_init_molecules_um_2'
	1.49E+02; %13.52143291;		% param(65) is 'kSon'
	7.50E-01; %4.6;		% param(66) is 'kG2SOSr'
	3.66E-01; %0.00942;		% param(67) is 'kG2SOSf'
	0.0;		% param(68) is 'Kr_neg_fdbk_dephos_iBRAF'
	0.0;		% param(69) is 'nfpRaf1_init_molecules_um_3'
	1.0;		% param(70) is 'netValence_Raf1_iBRaf_phos'
	1.43E-01; %0.072015388;		% param(71) is 'kiR1r'
	1.0;		% param(72) is 'netValence_EG2_SOS_bind'
	1.0;		% param(73) is 'netValence_neg_fdbk_unbind_RASt_nfpBRAF'
	300.0;		% param(74) is 'mlabfix_T_'
	1.58E-06; %1.15E-5;		% param(75) is 'kiR1f'
	1.0;		% param(76) is 'Size_cytoplasm'
	1.95E+02; %1243.069312;		% param(77) is 'kSoff'
	0.0;		% param(78) is 'Kr_MEK_dephos'
	0.0;		% param(79) is 'iBRaf_init_molecules_um_3'
	1.0;		% param(80) is 'netValence_BRaf_pRaf1_unbind'
	8314.46261815;		% param(81) is 'mlabfix_R_'
	0.0;		% param(82) is 'iBRaf_dimer_init_molecules_um_2'
	1.0;		% param(83) is 'netValence_BRAF_nfp'
	0.0;		% param(84) is 'mELmEL_init_molecules_um_2'
	0.0;		% param(85) is 'Ras_Braf_iRaf1_tetramer_init_molecules_um_2'
	1.33E-05; %6.0E-5;		% param(86) is 'kpERK'
	1.0;		% param(87) is 'netValence_Ras_iRaf1_bind'
	0.0;		% param(88) is 'Kr_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(89) is 'netValence_Ras_pRaf1_dp'
	0.0;		% param(90) is 'Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2'
	1.0;		% param(91) is 'netValence_EGF_bind'
	0.0;		% param(92) is 'Kr_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(93) is 'netValence_mELmEL_phos_dephos'
	0.0017;		% param(94) is 'EGF'
	2.69E-01; %0.142907544;		% param(95) is 'kRhydro'
	0.0;		% param(96) is 'Voltage_pm'
	0.0;		% param(97) is 'mEL_init_molecules_um_2'
	0.0;		% param(98) is 'Kr_pRAF1_RAF1_phos'
	0.0;		% param(99) is 'Kr_BRAF_nfp'
	1.0;		% param(100) is 'netValence_pRaf1_iBRaf_unbind'
	1.0;		% param(101) is 'netValence_Ras_Raf1_nfp'
	1.40E+05; %141781.6814;		% param(102) is 'Kmgneslow'
	4.31E-01; %0.6;		% param(103) is 'kdpMEK'
	2.08E-02; %0.046229858;		% param(104) is 'knfpiR1r'
	2.20E+01; %110.0;		% param(105) is 'kdR1r'
	2.51E-02; %0.0367;		% param(106) is 'kiBr'
	2.27E+00; %0.019;		% param(107) is 'kdR1f'
	7.68E-03; %2.47E-7;		% param(108) is 'kiBf'
	628000.0;		% param(109) is 'init_GRB2_HeLa'
	0.0;		% param(110) is 'Ras_nfpBRaf_init_molecules_um_2'
	2.94E-02; %0.37;		% param(111) is 'kfpBr'
	1.0;		% param(112) is 'netValence_BRaf_iBRaf_dim'
	0.0;		% param(113) is 'EG2SOS_init_molecules_um_2'
	4.02E-02; %0.6;		% param(114) is 'kdpSOS'
	12000.0;		% param(115) is 'init_RAF_HeLa'
	1.0;		% param(116) is 'netValence_iRaf_iBRaf_dim'
	0.0;		% param(117) is 'Kr_neg_fdbk_dephos_pSOS'
	1.0;		% param(118) is 'netValence_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(119) is 'netValence_iBRaf_Raf1_bind'
	3.72E+01; %8.0;		% param(120) is 'kdp'
	1.0;		% param(121) is 'netValence_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(122) is 'netValence_RAS_GTP_bind'
	0.0;		% param(123) is 'g0'
	1.0;		% param(124) is 'netValence_Ras_iRaf_nfp'
	0.0;		% param(125) is 'nfpBRaf_init_molecules_um_3'
	1.0;		% param(126) is 'netValence_Ras_Raf1_bind'
	1.0;		% param(127) is 'netValence_BRAF_RAF1_phos'
	1.61E-02; %0.1;		% param(128) is 'kdEr'
	3.141592653589793;		% param(129) is 'mlabfix_PI_'
	1.09E-01; %9.2E-4;		% param(130) is 'kdEf'
	3.11E+01; %0.6;		% param(131) is 'kdpR1'
	1000.0;		% param(132) is 'init_BRAF_HeLa'
	0.0;		% param(133) is 'BRaf_dimer_init_molecules_um_2'
	0.0;		% param(134) is 'Kr_neg_fdbck_dephos_piRaf'
	0.0;		% param(135) is 'Kr_neg_fdbk_phos_of_SOS'
	0.0;		% param(136) is 'Kr_neg_fdbk_unbind_RASt_nfpRAF1'
	0.0;		% param(137) is 'Ras_Raf1_iBRaf_tetramer_init_molecules_um_2'
	0.0;		% param(138) is 'BRaf_iBRaf_dimer_init_molecules_um_2'
	135000.0;		% param(139) is 'init_RAS_HeLa'
	6.02214179E11;		% param(140) is 'mlabfix_N_pmol_'
	0.0;		% param(141) is 'Ras_pRaf1_init_molecules_um_2'
	93000.0;		% param(142) is 'init_EGFR_HeLa'
	0.0;		% param(143) is 'nfpiRaf1_init_molecules_um_3'
	2.18E+00; %460.0;		% param(144) is 'kG2r'
	1.0;		% param(145) is 'netValence_RAS_GTPhydro'
	0.0;		% param(146) is 'Ras_nfpiRaf1_init_molecules_um_2'
	10.0;		% param(147) is 'init_sorafenib'
	4.66E-05; %3.8E-4;		% param(148) is 'kG2f'
	0.0;		% param(149) is 'Ras_Raf1_init_molecules_um_2'
	0.0;		% param(150) is 'Kr_Ras_pRaf1_dp'
	0.001660538783162726;		% param(151) is 'KMOLE'
	0.0;		% param(152) is 'Kr_Raf1_iBRaf_phos'
	0.0;		% param(153) is 'Kr_neg_fdbk_dephos_pBRAF'
	0.0;		% param(154) is 'pERK_init_molecules_um_3'
	0.0;		% param(155) is 'Ras_nfpiBRaf_init_molecules_um_2'
	1.0;		% param(156) is 'Size_pm'
	1.0;		% param(157) is 'netValence_BRAF_dim'
	0.0;		% param(158) is 'nfpSOS_init_molecules_um_3'
	0.0;		% param(159) is 'Ras_BRaf_Raf1_tetramer_init_molecules_um_2'
	1.60E+00; %0.37;		% param(160) is 'knfpiBr'
	0.0;		% param(161) is 'Kr_Ras_Raf1_nfp'
	0.0;		% param(162) is 'Ras_iRaf1_init_molecules_um_2'
	7.80E+02; %0.046229858;		% param(163) is 'kfpR1r'
	3.80E-02; %0.8;		% param(164) is 'kEr'
];
%}


if strcmp(modeltype, 'min_sor') == 1
	param(146) = 0;
    
    %param(65) = 0; %kSon
    %param(71) = 0; % kiR1r
    %param(75) = 0; % kiR1f
    %param(106) = 0; %kiBr
    %param(108) = 0; %kiBf
    %param(104) = 0; %knfpiR1r
    %param(159) = 0; %knfpiBr
end

if strcmp(modeltype, 'plus_sor') == 1
	param(146) = 10;
end
%
% invoke the integrator
%[T,Y] = ode15s(@f,timeSpan,yinit,odeset,param,yinit,options);
%[T,Y] = ode15s(@f,timeSpan,yinit,odeset('RelTol',1e-9,'AbsTol',1e-9),param,yinit);
options = odeset('RelTol',1e-9,'AbsTol',1e-9);
%[T,Y] = ode15s(@f,timeSpan,yinit,odeset('RelTol',1e-9,'AbsTol',1e-9),param);
[T,Y] = ode15s(@f,timeSpan,yinit,odeset('RelTol',1e-9,'AbsTol',1e-9),param,yinit);


%f(t,y,p,y0)
%[T,Y] = ode15s(@f, timeSpan, yinit, options, param, yinit);


% get the solution
all = zeros(length(T), numVars);
for i = 1:length(T)
	all(i,:) = getRow(T(i), Y(i,:), yinit, param);
end

allValues = all;
end

% -------------------------------------------------------
% get row data
function rowValue = getRow(t,y,y0,p)
	% State Variables
	Ras_iBRaf = y(1);
	Ras_Raf1_iBRaf_tetramer = y(2);
	Ras_BRaf = y(3);
	Ras_nfpRaf1 = y(4);
	iRaf1 = y(5);
	ERK = y(6);
	nfpSOS = y(7);
	Ras_pRaf1 = y(8);
	BRaf_iBRaf_dimer = y(9);
	Ras_pRaf1_iBRaf_tetramer = y(10);
	Ras_pRaf1_tetramer = y(11);
	Ras_nfpiRaf1 = y(12);
	Ras_GTP = y(13);
	nfpiRaf1 = y(14);
	Raf1 = y(15);
	Ras_BRaf_pRaf1_tetramer = y(16);
	Ras_pRaf1_Raf1_tetramer = y(17);
	nfpRaf1 = y(18);
	Ras_nfpBRaf = y(19);
	iBRaf = y(20);
	pMEK = y(21);
	mE = y(22);
	mEL = y(23);
	pRaf1 = y(24);
	EG2 = y(25);
	Ras_iRaf1 = y(26);
	Ras_nfpiBRaf = y(27);
	MEK = y(28);
	Ras_Braf_iRaf1_tetramer = y(29);
	Ras_Raf1 = y(30);
	mELmEL = y(31);
	nfpiBRaf = y(32);
	BRaf = y(33);
	SOS = y(34);
	nfpBRaf = y(35);
	GRB2_SOS = y(36);
	GRB2 = y(37);
	Ras_iRaf1_tetramer = y(38);
	iBRaf_dimer = y(39);
	Ras_GDP = y(40);
	E = y(41);
	Ras_BRaf_Raf1_tetramer = y(42);
	BRaf_dimer = y(43);
	EG2SOS = y(44);
	pERK = y(45);
	Ras_iRaf1_iBRaf1_tetramer = y(46);
	% Constants
	Ras_iRaf1_tetramer_init_molecules_um_2 = p(1);
	kEf = p(2);
	mlabfix_F_ = p(3);
	Ras_pRaf1_Raf1_tetramer_init_molecules_um_2 = p(4);
	kcatE = p(5);
	pRaf1_init_molecules_um_3 = p(6);
	netValence_neg_fdbk_phos_of_SOS = p(7);
	kpMEK = p(8);
	netValence_E_G2SOS_bind = p(9);
	knfpSOS = p(10);
	init_SOS_HeLa = p(11);
	netValence_pRaf1_Raf1_bind = p(12);
	netValence_iBRaf_dim = p(13);
	netValence_BRaf_iRaf1_bind = p(14);
	netValence_neg_fdbk_unbind_RASt_nfpRAF1 = p(15);
	kBr = p(16);
	Kr_ERK_phos = p(17);
	Kr_RAS_GTP_bind = p(18);
	init_ERK_HeLa = p(19);
	Kr_Ras_iRaf_nfp = p(20);
	netValence_pRAF1_RAF1_phos = p(21);
	kBf = p(22);
	kdpERK = p(23);
	Kr_pRAF1_dephos = p(24);
	kR1r = p(25);
	iRaf1_init_molecules_um_3 = p(26);
	kR1f = p(27);
	netValence_RAS_pRAF1_unbind = p(28);
	Kr_neg_fdbk_unbind_RASt_nfpBRAF = p(29);
	Kr_MEK_phos = p(30);
	init_MEK_HeLa = p(31);
	netValence_Ras_iBRaf_nfp = p(32);
	Kr_neg_fdbk_dephos_pRaf1 = p(33);
	netValence_RAS_iBRAF_bind = p(34);
	nfpiBRaf_init_molecules_um_3 = p(35);
	netValence_E_G2_bind = p(36);
	Kr_RAS_GTPhydro = p(37);
	Kr_Ras_iBRaf_nfp = p(38);
	pMEK_init_molecules_um_3 = p(39);
	Ras_BRaf_init_molecules_um_2 = p(40);
	Ras_iRaf1_iBRaf1_tetramer_init_molecules_um_2 = p(41);
	Ras_pRaf1_tetramer_init_molecules_um_2 = p(42);
	K_millivolts_per_volt = p(43);
	mlabfix_K_GHK_ = p(44);
	knfpBR1 = p(45);
	EG2_init_molecules_um_2 = p(46);
	Ras_BRaf_pRaf1_tetramer_init_molecules_um_2 = p(47);
	netValence_pRaf1_pRaf1_unbind = p(48);
	kpR1 = p(49);
	kSOSr = p(50);
	Ras_nfpRaf1_init_molecules_um_2 = p(51);
	netValence_iRaf1_iRaf1_bind = p(52);
	kSOSf = p(53);
	Ras_iBRaf_init_molecules_um_2 = p(54);
	netValence_RAS_BRAF_bind = p(55);
	Kr_ERK_dephos = p(56);
	kRgneslow = p(57);
	mlabfix_F_nmol_ = p(58);
	netValence_mEL_dim = p(59);
	netValence_BRAF_RAF1_bind = p(60);
	GRB2_SOS_init_molecules_um_3 = p(61);
	Kr_BRAF_RAF1_phos = p(62);
	E_init_molecules_um_2 = p(63);
	Ras_GTP_init_molecules_um_2 = p(64);
	kSon = p(65);
	kG2SOSr = p(66);
	kG2SOSf = p(67);
	Kr_neg_fdbk_dephos_iBRAF = p(68);
	nfpRaf1_init_molecules_um_3 = p(69);
	netValence_Raf1_iBRaf_phos = p(70);
	kiR1r = p(71);
	netValence_EG2_SOS_bind = p(72);
	netValence_neg_fdbk_unbind_RASt_nfpBRAF = p(73);
	mlabfix_T_ = p(74);
	kiR1f = p(75);
	Size_cytoplasm = p(76);
	kSoff = p(77);
	Kr_MEK_dephos = p(78);
	iBRaf_init_molecules_um_3 = p(79);
	netValence_BRaf_pRaf1_unbind = p(80);
	mlabfix_R_ = p(81);
	iBRaf_dimer_init_molecules_um_2 = p(82);
	netValence_BRAF_nfp = p(83);
	mELmEL_init_molecules_um_2 = p(84);
	Ras_Braf_iRaf1_tetramer_init_molecules_um_2 = p(85);
	kpERK = p(86);
	netValence_Ras_iRaf1_bind = p(87);
	Kr_neg_fdbk_unbind_piRaf1_RASt = p(88);
	netValence_Ras_pRaf1_dp = p(89);
	Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2 = p(90);
	netValence_EGF_bind = p(91);
	Kr_neg_fdbk_unbind_piBRaf_RASt = p(92);
	netValence_mELmEL_phos_dephos = p(93);
	EGF = p(94);
	kRhydro = p(95);
	Voltage_pm = p(96);
	mEL_init_molecules_um_2 = p(97);
	Kr_pRAF1_RAF1_phos = p(98);
	Kr_BRAF_nfp = p(99);
	netValence_pRaf1_iBRaf_unbind = p(100);
	netValence_Ras_Raf1_nfp = p(101);
	Kmgneslow = p(102);
	kdpMEK = p(103);
	knfpiR1r = p(104);
	kdR1r = p(105);
	kiBr = p(106);
	kdR1f = p(107);
	kiBf = p(108);
	init_GRB2_HeLa = p(109);
	Ras_nfpBRaf_init_molecules_um_2 = p(110);
	kfpBr = p(111);
	netValence_BRaf_iBRaf_dim = p(112);
	EG2SOS_init_molecules_um_2 = p(113);
	kdpSOS = p(114);
	init_RAF_HeLa = p(115);
	netValence_iRaf_iBRaf_dim = p(116);
	Kr_neg_fdbk_dephos_pSOS = p(117);
	netValence_neg_fdbk_unbind_piRaf1_RASt = p(118);
	netValence_iBRaf_Raf1_bind = p(119);
	kdp = p(120);
	netValence_neg_fdbk_unbind_piBRaf_RASt = p(121);
	netValence_RAS_GTP_bind = p(122);
	netValence_Ras_iRaf_nfp = p(123);
	nfpBRaf_init_molecules_um_3 = p(124);
	netValence_Ras_Raf1_bind = p(125);
	netValence_BRAF_RAF1_phos = p(126);
	kdEr = p(127);
	mlabfix_PI_ = p(128);
	kdEf = p(129);
	kdpR1 = p(130);
	init_BRAF_HeLa = p(131);
	BRaf_dimer_init_molecules_um_2 = p(132);
	Kr_neg_fdbck_dephos_piRaf = p(133);
	Kr_neg_fdbk_phos_of_SOS = p(134);
	Kr_neg_fdbk_unbind_RASt_nfpRAF1 = p(135);
	Ras_Raf1_iBRaf_tetramer_init_molecules_um_2 = p(136);
	BRaf_iBRaf_dimer_init_molecules_um_2 = p(137);
	init_RAS_HeLa = p(138);
	mlabfix_N_pmol_ = p(139);
	Ras_pRaf1_init_molecules_um_2 = p(140);
	init_EGFR_HeLa = p(141);
	nfpiRaf1_init_molecules_um_3 = p(142);
	kG2r = p(143);
	netValence_RAS_GTPhydro = p(144);
	Ras_nfpiRaf1_init_molecules_um_2 = p(145);
	init_sorafenib = p(146);
	kG2f = p(147);
	Ras_Raf1_init_molecules_um_2 = p(148);
	Kr_Ras_pRaf1_dp = p(149);
	KMOLE = p(150);
	Kr_Raf1_iBRaf_phos = p(151);
	Kr_neg_fdbk_dephos_pBRAF = p(152);
	pERK_init_molecules_um_3 = p(153);
	Ras_nfpiBRaf_init_molecules_um_2 = p(154);
	Size_pm = p(155);
	netValence_BRAF_dim = p(156);
	nfpSOS_init_molecules_um_3 = p(157);
	Ras_BRaf_Raf1_tetramer_init_molecules_um_2 = p(158);
	knfpiBr = p(159);
	Kr_Ras_Raf1_nfp = p(160);
	Ras_iRaf1_init_molecules_um_2 = p(161);
	kfpR1r = p(162);
	kEr = p(163);
	% Functions
	Kr_E_G2_bind = kG2r;
	Ras_GDP_init_molecules_um_2 = init_RAS_HeLa;
	Kf_iRaf_iBRaf_dim = kdR1f;
	Raf1_cyt_param = (Raf1 + pRaf1 + iRaf1 + nfpRaf1 + nfpiRaf1);
	Kr_E_G2SOS_bind = kG2SOSr;
	Kf_E_G2SOS_bind = kG2SOSf;
	J_E_G2SOS_bind = (((Kf_E_G2SOS_bind .* GRB2_SOS) .* E) - (Kr_E_G2SOS_bind .* EG2SOS));
	Kr_BRAF_dim = kdR1r;
	Kf_RAS_GTP_bind = (kRgneslow .* EG2SOS ./ (Kmgneslow + Ras_GDP));
	J_RAS_GTP_bind = ((Kf_RAS_GTP_bind .* Ras_GDP) - (Kr_RAS_GTP_bind .* Ras_GTP));
	Kf_Ras_iRaf_nfp = (kpR1 .* pERK);
	J_Ras_iRaf_nfp = ((Kf_Ras_iRaf_nfp .* Ras_iRaf1) - (Kr_Ras_iRaf_nfp .* Ras_nfpiRaf1));
	Kf_Ras_Raf1_nfp = (knfpBR1 .* pERK);
	Kf_Ras_iBRaf_nfp = (knfpBR1 .* pERK);
	Kf_pRAF1_dephos = kdpR1;
	J_pRAF1_dephos = ((Kf_pRAF1_dephos .* pRaf1) - (Kr_pRAF1_dephos .* Raf1));
	Kf_iRaf1_iRaf1_bind = kdR1f;
	Kf_BRAF_RAF1_phos = kpR1;
	Kf_BRAF_dim = kdR1f;
	J_BRAF_dim = ((Kf_BRAF_dim .* (Ras_BRaf ^ 2.0)) - (Kr_BRAF_dim .* BRaf_dimer));
	sorafenib_init_molecules_um_3 = init_sorafenib;
	sorafenib = sorafenib_init_molecules_um_3;
	Kf_E_G2_bind = kG2f;
	J_E_G2_bind = (((Kf_E_G2_bind .* GRB2) .* E) - (Kr_E_G2_bind .* EG2));
	Kf_BRAF_nfp = (knfpBR1 .* pERK);
	Kf_neg_fdbk_dephos_pSOS = kdpSOS;
	J_neg_fdbk_dephos_pSOS = ((Kf_neg_fdbk_dephos_pSOS .* nfpSOS) - (Kr_neg_fdbk_dephos_pSOS .* SOS));
	Kr_Ras_Raf1_bind = kR1r;
	Kf_neg_fdbk_dephos_iBRAF = kdpR1;
	Kf_ERK_dephos = kdpERK;
	Kf_iBRaf_dim = kdR1f;
	pRaf1_total = (Ras_pRaf1 + Ras_pRaf1_iBRaf_tetramer + Ras_pRaf1_tetramer + Ras_BRaf_pRaf1_tetramer + pRaf1 + Ras_pRaf1_Raf1_tetramer);
	Kf_EGF_bind = (kEf .* EGF);
	Kf_EG2_SOS_bind = kSOSf;
	Kf_RAS_GTPhydro = kRhydro;
	J_RAS_GTPhydro = ((Kf_RAS_GTPhydro .* Ras_GTP) - (Kr_RAS_GTPhydro .* Ras_GDP));
	Kf_iBRaf_Raf1_bind = kdR1f;
	Kf_neg_fdbk_unbind_RASt_nfpBRAF = kfpBr;
	Kr_mELmEL_phos_dephos = kdp;
	Kf_RAS_BRAF_bind = kBf;
	Kr_iBRaf_Raf1_bind = kdR1r;
	Kf_pRaf1_Raf1_bind = kdR1f;
	Kf_neg_fdbk_phos_of_SOS = (knfpSOS .* pERK);
	Kr_RAS_pRAF1_unbind = kR1f;
	Kf_neg_fdbk_unbind_piRaf1_RASt = knfpiR1r;
	Kf_Ras_Raf1_bind = kR1f;
	J_Ras_Raf1_bind = (((Kf_Ras_Raf1_bind .* Ras_GTP) .* Raf1) - (Kr_Ras_Raf1_bind .* Ras_Raf1));
	Kf_neg_fdbk_unbind_piBRaf_RASt = knfpiBr;
	Kf_Raf1_iBRaf_phos = kpR1;
	J_Raf1_iBRaf_phos = ((Kf_Raf1_iBRaf_phos .* Ras_Raf1_iBRaf_tetramer) - (Kr_Raf1_iBRaf_phos .* Ras_pRaf1_iBRaf_tetramer));
	Kr_BRaf_iBRaf_dim = kdR1r;
	Kf_BRaf_iBRaf_dim = kdR1f;
	Kf_sorafenib_Raf1_bind_cyt = kSon;
	J_ERK_dephos = ((Kf_ERK_dephos .* pERK) - (Kr_ERK_dephos .* ERK));
	Kf_Ras_pRaf1_dp = kdpR1;
	Kf_GRB2_SOS_bind_cyt = kSOSf;
	Kr_GRB2_SOS_bind_cyt = kSOSr;
	J_GRB2_SOS_bind_cyt = (((Kf_GRB2_SOS_bind_cyt .* GRB2) .* SOS) - (Kr_GRB2_SOS_bind_cyt .* GRB2_SOS));
	Kf_mEL_dim = kdEf;
	Kf_MEK_dephos = kdpMEK;
	J_neg_fdbk_unbind_RASt_nfpBRAF = ((Kf_neg_fdbk_unbind_RASt_nfpBRAF .* Ras_nfpBRaf) - ((Kr_neg_fdbk_unbind_RASt_nfpBRAF .* Ras_GTP) .* nfpBRaf));
	Kr_iRaf_iBRaf_dim = kdR1r;
	Kr_RAS_iBRAF_bind = kiBr;
	Kf_RAS_iBRAF_bind = kiBf;
	J_RAS_iBRAF_bind = (((Kf_RAS_iBRAF_bind .* iBRaf) .* Ras_GTP) - (Kr_RAS_iBRAF_bind .* Ras_iBRaf));
	Kf_neg_fdbk_dephos_pBRAF = kdpR1;
	J_neg_fdbk_dephos_pBRAF = ((Kf_neg_fdbk_dephos_pBRAF .* nfpBRaf) - (Kr_neg_fdbk_dephos_pBRAF .* BRaf));
	Kf_Ras_iRaf1_bind = kiR1f;
	Kr_Ras_iRaf1_bind = kiR1r;
	J_Ras_iRaf1_bind = (((Kf_Ras_iRaf1_bind .* iRaf1) .* Ras_GTP) - (Kr_Ras_iRaf1_bind .* Ras_iRaf1));
	SOS_init_molecules_um_3 = init_SOS_HeLa;
	mE_init_molecules_um_2 = init_EGFR_HeLa;
	KFlux_pm_cytoplasm = (Size_pm ./ Size_cytoplasm);
	ERK_init_molecules_um_3 = init_ERK_HeLa;
	Kr_EG2_SOS_bind = kSOSr;
	J_EG2_SOS_bind = (((Kf_EG2_SOS_bind .* EG2) .* SOS) - (Kr_EG2_SOS_bind .* EG2SOS));
	J_MEK_dephos = ((Kf_MEK_dephos .* pMEK) - (Kr_MEK_dephos .* MEK));
	Kf_ERK_phos = (kpERK .* pMEK);
	J_ERK_phos = ((Kf_ERK_phos .* ERK) - (Kr_ERK_phos .* pERK));
	Kf_BRAF_RAF1_bind = kdR1f;
	Kr_BRAF_RAF1_bind = kdR1r;
	J_BRAF_RAF1_bind = (((Kf_BRAF_RAF1_bind .* Ras_BRaf) .* Ras_Raf1) - (Kr_BRAF_RAF1_bind .* Ras_BRaf_Raf1_tetramer));
	Kf_pRaf1_pRaf1_unbind = kdR1r;
	Kr_pRaf1_pRaf1_unbind = kdR1f;
	J_pRaf1_pRaf1_unbind = ((Kf_pRaf1_pRaf1_unbind .* Ras_pRaf1_tetramer) - (Kr_pRaf1_pRaf1_unbind .* (Ras_pRaf1 ^ 2.0)));
	active_Raf = ((2.0 .* BRaf_dimer) + Ras_BRaf_Raf1_tetramer + (2.0 .* Ras_BRaf_pRaf1_tetramer) + Ras_pRaf1_iBRaf_tetramer + Ras_pRaf1 + Ras_pRaf1_Raf1_tetramer + (2.0 .* Ras_pRaf1_tetramer) + pRaf1 + BRaf_iBRaf_dimer + Ras_Braf_iRaf1_tetramer);
	Kf_MEK_phos = (kpMEK .* active_Raf);
	J_MEK_phos = ((Kf_MEK_phos .* MEK) - (Kr_MEK_phos .* pMEK));
	Kf_BRaf_pRaf1_unbind = kdR1r;
	Kr_BRaf_pRaf1_unbind = kdR1f;
	J_BRaf_pRaf1_unbind = ((Kf_BRaf_pRaf1_unbind .* Ras_BRaf_pRaf1_tetramer) - ((Kr_BRaf_pRaf1_unbind .* Ras_BRaf) .* Ras_pRaf1));
	Kf_BRaf_iRaf1_bind = kdR1f;
	Kf_pRaf1_iBRaf_unbind = kdR1r;
	Kr_pRaf1_iBRaf_unbind = kdR1f;
	J_pRaf1_iBRaf_unbind = ((Kf_pRaf1_iBRaf_unbind .* Ras_pRaf1_iBRaf_tetramer) - ((Kr_pRaf1_iBRaf_unbind .* Ras_iBRaf) .* Ras_pRaf1));
	Kf_sorafenib_BRaf_bind_cyt = kSon;
	J_iBRaf_Raf1_bind = (((Kf_iBRaf_Raf1_bind .* Ras_Raf1) .* Ras_iBRaf) - (Kr_iBRaf_Raf1_bind .* Ras_Raf1_iBRaf_tetramer));
	Kr_mEL_dim = kdEr;
	J_mEL_dim = ((Kf_mEL_dim .* (mEL ^ 2.0)) - (Kr_mEL_dim .* mELmEL));
	Kr_pRaf1_Raf1_bind = kdR1r;
	Kf_pRAF1_RAF1_phos = kpR1;
	Kf_RAS_pRAF1_unbind = kR1r;
	J_RAS_pRAF1_unbind = ((Kf_RAS_pRAF1_unbind .* Ras_pRaf1) - ((Kr_RAS_pRAF1_unbind .* pRaf1) .* Ras_GTP));
	Kr_BRaf_iRaf1_bind = kdR1r;
	Kf_neg_fdbck_dephos_piRaf = kdpR1;
	J_Ras_Raf1_nfp = ((Kf_Ras_Raf1_nfp .* Ras_Raf1) - (Kr_Ras_Raf1_nfp .* Ras_nfpRaf1));
	Kf_neg_fdbk_dephos_pRaf1 = kdpR1;
	J_neg_fdbk_dephos_pRaf1 = ((Kf_neg_fdbk_dephos_pRaf1 .* nfpRaf1) - (Kr_neg_fdbk_dephos_pRaf1 .* Raf1));
	Kr_sorafenib_Raf1_bind_cyt = kSoff;
	J_sorafenib_Raf1_bind_cyt = (((Kf_sorafenib_Raf1_bind_cyt .* sorafenib) .* Raf1) - (Kr_sorafenib_Raf1_bind_cyt .* iRaf1));
	Raf1_pm_param = (Ras_pRaf1 + (2.0 .* Ras_pRaf1_Raf1_tetramer) + Ras_BRaf_Raf1_tetramer + Ras_iRaf1 + Ras_Raf1 + (2.0 .* Ras_iRaf1_tetramer) + Ras_Braf_iRaf1_tetramer + Ras_BRaf_pRaf1_tetramer + (2.0 .* Ras_pRaf1_tetramer) + Ras_nfpiRaf1 + Ras_nfpRaf1 + Ras_iRaf1_iBRaf1_tetramer + Ras_Raf1_iBRaf_tetramer + Ras_pRaf1_iBRaf_tetramer);
	boundRaf_frac_parameter = (Raf1_pm_param ./ (Raf1_cyt_param + Raf1_pm_param));
	Kr_EGF_bind = kEr;
	J_BRAF_nfp = ((Kf_BRAF_nfp .* Ras_BRaf) - (Kr_BRAF_nfp .* Ras_nfpBRaf));
	J_Ras_iBRaf_nfp = ((Kf_Ras_iBRaf_nfp .* Ras_iBRaf) - (Kr_Ras_iBRaf_nfp .* Ras_nfpiBRaf));
	Kr_RAS_BRAF_bind = kBr;
	J_EGF_bind = ((Kf_EGF_bind .* mE) - (Kr_EGF_bind .* mEL));
	Kr_iRaf1_iRaf1_bind = kdR1r;
	MEK_init_molecules_um_3 = init_MEK_HeLa;
	J_iRaf1_iRaf1_bind = ((Kf_iRaf1_iRaf1_bind .* (Ras_iRaf1 ^ 2.0)) - (Kr_iRaf1_iRaf1_bind .* Ras_iRaf1_tetramer));
	GRB2_init_molecules_um_3 = init_GRB2_HeLa;
	J_BRaf_iBRaf_dim = (((Kf_BRaf_iBRaf_dim .* Ras_iBRaf) .* Ras_BRaf) - (Kr_BRaf_iBRaf_dim .* BRaf_iBRaf_dimer));
	Kf_mELmEL_phos_dephos = kcatE;
	J_mELmEL_phos_dephos = ((Kf_mELmEL_phos_dephos .* mELmEL) - (Kr_mELmEL_phos_dephos .* E));
	J_RAS_BRAF_bind = (((Kf_RAS_BRAF_bind .* BRaf) .* Ras_GTP) - (Kr_RAS_BRAF_bind .* Ras_BRaf));
	J_iRaf_iBRaf_dim = (((Kf_iRaf_iBRaf_dim .* Ras_iBRaf) .* Ras_iRaf1) - (Kr_iRaf_iBRaf_dim .* Ras_iRaf1_iBRaf1_tetramer));
	J_neg_fdbk_dephos_iBRAF = ((Kf_neg_fdbk_dephos_iBRAF .* nfpiBRaf) - (Kr_neg_fdbk_dephos_iBRAF .* iBRaf));
	J_Ras_pRaf1_dp = ((Kf_Ras_pRaf1_dp .* Ras_pRaf1) - (Kr_Ras_pRaf1_dp .* Ras_Raf1));
	Kr_iBRaf_dim = kdR1r;
	Kf_neg_fdbk_unbind_RASt_nfpRAF1 = kfpR1r;
	Raf1_init_molecules_um_3 = init_RAF_HeLa;
	J_BRAF_RAF1_phos = ((Kf_BRAF_RAF1_phos .* Ras_BRaf_Raf1_tetramer) - (Kr_BRAF_RAF1_phos .* Ras_BRaf_pRaf1_tetramer));
	Kr_sorafenib_BRaf_bind_cyt = kSoff;
	J_sorafenib_BRaf_bind_cyt = (((Kf_sorafenib_BRaf_bind_cyt .* sorafenib) .* BRaf) - (Kr_sorafenib_BRaf_bind_cyt .* iBRaf));
	J_pRaf1_Raf1_bind = (((Kf_pRaf1_Raf1_bind .* Ras_pRaf1) .* Ras_Raf1) - (Kr_pRaf1_Raf1_bind .* Ras_pRaf1_Raf1_tetramer));
	J_neg_fdbk_unbind_piRaf1_RASt = ((Kf_neg_fdbk_unbind_piRaf1_RASt .* Ras_nfpiRaf1) - ((Kr_neg_fdbk_unbind_piRaf1_RASt .* nfpiRaf1) .* Ras_GTP));
	BRaf_init_molecules_um_3 = init_BRAF_HeLa;
	RasGTP_parameter = (Ras_GTP + (2.0 .* (Ras_iRaf1_tetramer + Ras_Braf_iRaf1_tetramer + Ras_BRaf_Raf1_tetramer + Ras_pRaf1_Raf1_tetramer + Ras_BRaf_pRaf1_tetramer + Ras_pRaf1_tetramer + BRaf_dimer + BRaf_iBRaf_dimer + iBRaf_dimer + Ras_iRaf1_iBRaf1_tetramer + Ras_Raf1_iBRaf_tetramer + Ras_pRaf1_iBRaf_tetramer)) + Ras_BRaf + Ras_iRaf1 + Ras_pRaf1 + Ras_Raf1 + Ras_nfpRaf1 + Ras_nfpiRaf1 + Ras_iBRaf + Ras_nfpiBRaf + Ras_nfpBRaf);
	J_BRaf_iRaf1_bind = (((Kf_BRaf_iRaf1_bind .* Ras_BRaf) .* Ras_iRaf1) - (Kr_BRaf_iRaf1_bind .* Ras_Braf_iRaf1_tetramer));
	J_neg_fdbk_phos_of_SOS = ((Kf_neg_fdbk_phos_of_SOS .* EG2SOS) - ((Kr_neg_fdbk_phos_of_SOS .* EG2) .* nfpSOS));
	J_neg_fdbk_unbind_piBRaf_RASt = ((Kf_neg_fdbk_unbind_piBRaf_RASt .* Ras_nfpiBRaf) - ((Kr_neg_fdbk_unbind_piBRaf_RASt .* nfpiBRaf) .* Ras_GTP));
	J_iBRaf_dim = ((Kf_iBRaf_dim .* (Ras_iBRaf ^ 2.0)) - (Kr_iBRaf_dim .* iBRaf_dimer));
	J_neg_fdbck_dephos_piRaf = ((Kf_neg_fdbck_dephos_piRaf .* nfpiRaf1) - (Kr_neg_fdbck_dephos_piRaf .* iRaf1));
	J_pRAF1_RAF1_phos = ((Kf_pRAF1_RAF1_phos .* Ras_pRaf1_Raf1_tetramer) - (Kr_pRAF1_RAF1_phos .* Ras_pRaf1_tetramer));
	J_neg_fdbk_unbind_RASt_nfpRAF1 = ((Kf_neg_fdbk_unbind_RASt_nfpRAF1 .* Ras_nfpRaf1) - ((Kr_neg_fdbk_unbind_RASt_nfpRAF1 .* Ras_GTP) .* nfpRaf1));
	% OutputFunctions

	rowValue = [Ras_iBRaf Ras_Raf1_iBRaf_tetramer Ras_BRaf Ras_nfpRaf1 iRaf1 ERK nfpSOS Ras_pRaf1 BRaf_iBRaf_dimer Ras_pRaf1_iBRaf_tetramer Ras_pRaf1_tetramer Ras_nfpiRaf1 Ras_GTP nfpiRaf1 Raf1 Ras_BRaf_pRaf1_tetramer Ras_pRaf1_Raf1_tetramer nfpRaf1 Ras_nfpBRaf iBRaf pMEK mE mEL pRaf1 EG2 Ras_iRaf1 Ras_nfpiBRaf MEK Ras_Braf_iRaf1_tetramer Ras_Raf1 mELmEL nfpiBRaf BRaf SOS nfpBRaf GRB2_SOS GRB2 Ras_iRaf1_tetramer iBRaf_dimer Ras_GDP E Ras_BRaf_Raf1_tetramer BRaf_dimer EG2SOS pERK Ras_iRaf1_iBRaf1_tetramer Kr_E_G2_bind Ras_GDP_init_molecules_um_2 Kf_iRaf_iBRaf_dim Raf1_cyt_param Kr_E_G2SOS_bind Kf_E_G2SOS_bind J_E_G2SOS_bind Kr_BRAF_dim Kf_RAS_GTP_bind J_RAS_GTP_bind Kf_Ras_iRaf_nfp J_Ras_iRaf_nfp Kf_Ras_Raf1_nfp Kf_Ras_iBRaf_nfp Kf_pRAF1_dephos J_pRAF1_dephos Kf_iRaf1_iRaf1_bind Kf_BRAF_RAF1_phos Kf_BRAF_dim J_BRAF_dim sorafenib_init_molecules_um_3 sorafenib Kf_E_G2_bind J_E_G2_bind Kf_BRAF_nfp Kf_neg_fdbk_dephos_pSOS J_neg_fdbk_dephos_pSOS Kr_Ras_Raf1_bind Kf_neg_fdbk_dephos_iBRAF Kf_ERK_dephos Kf_iBRaf_dim pRaf1_total Kf_EGF_bind Kf_EG2_SOS_bind Kf_RAS_GTPhydro J_RAS_GTPhydro Kf_iBRaf_Raf1_bind Kf_neg_fdbk_unbind_RASt_nfpBRAF Kr_mELmEL_phos_dephos Kf_RAS_BRAF_bind Kr_iBRaf_Raf1_bind Kf_pRaf1_Raf1_bind Kf_neg_fdbk_phos_of_SOS Kr_RAS_pRAF1_unbind Kf_neg_fdbk_unbind_piRaf1_RASt Kf_Ras_Raf1_bind J_Ras_Raf1_bind Kf_neg_fdbk_unbind_piBRaf_RASt Kf_Raf1_iBRaf_phos J_Raf1_iBRaf_phos Kr_BRaf_iBRaf_dim Kf_BRaf_iBRaf_dim Kf_sorafenib_Raf1_bind_cyt J_ERK_dephos Kf_Ras_pRaf1_dp Kf_GRB2_SOS_bind_cyt Kr_GRB2_SOS_bind_cyt J_GRB2_SOS_bind_cyt Kf_mEL_dim Kf_MEK_dephos J_neg_fdbk_unbind_RASt_nfpBRAF Kr_iRaf_iBRaf_dim Kr_RAS_iBRAF_bind Kf_RAS_iBRAF_bind J_RAS_iBRAF_bind Kf_neg_fdbk_dephos_pBRAF J_neg_fdbk_dephos_pBRAF Kf_Ras_iRaf1_bind Kr_Ras_iRaf1_bind J_Ras_iRaf1_bind SOS_init_molecules_um_3 mE_init_molecules_um_2 KFlux_pm_cytoplasm ERK_init_molecules_um_3 Kr_EG2_SOS_bind J_EG2_SOS_bind J_MEK_dephos Kf_ERK_phos J_ERK_phos Kf_BRAF_RAF1_bind Kr_BRAF_RAF1_bind J_BRAF_RAF1_bind Kf_pRaf1_pRaf1_unbind Kr_pRaf1_pRaf1_unbind J_pRaf1_pRaf1_unbind active_Raf Kf_MEK_phos J_MEK_phos Kf_BRaf_pRaf1_unbind Kr_BRaf_pRaf1_unbind J_BRaf_pRaf1_unbind Kf_BRaf_iRaf1_bind Kf_pRaf1_iBRaf_unbind Kr_pRaf1_iBRaf_unbind J_pRaf1_iBRaf_unbind Kf_sorafenib_BRaf_bind_cyt J_iBRaf_Raf1_bind Kr_mEL_dim J_mEL_dim Kr_pRaf1_Raf1_bind Kf_pRAF1_RAF1_phos Kf_RAS_pRAF1_unbind J_RAS_pRAF1_unbind Kr_BRaf_iRaf1_bind Kf_neg_fdbck_dephos_piRaf J_Ras_Raf1_nfp Kf_neg_fdbk_dephos_pRaf1 J_neg_fdbk_dephos_pRaf1 Kr_sorafenib_Raf1_bind_cyt J_sorafenib_Raf1_bind_cyt Raf1_pm_param boundRaf_frac_parameter Kr_EGF_bind J_BRAF_nfp J_Ras_iBRaf_nfp Kr_RAS_BRAF_bind J_EGF_bind Kr_iRaf1_iRaf1_bind MEK_init_molecules_um_3 J_iRaf1_iRaf1_bind GRB2_init_molecules_um_3 J_BRaf_iBRaf_dim Kf_mELmEL_phos_dephos J_mELmEL_phos_dephos J_RAS_BRAF_bind J_iRaf_iBRaf_dim J_neg_fdbk_dephos_iBRAF J_Ras_pRaf1_dp Kr_iBRaf_dim Kf_neg_fdbk_unbind_RASt_nfpRAF1 Raf1_init_molecules_um_3 J_BRAF_RAF1_phos Kr_sorafenib_BRaf_bind_cyt J_sorafenib_BRaf_bind_cyt J_pRaf1_Raf1_bind J_neg_fdbk_unbind_piRaf1_RASt BRaf_init_molecules_um_3 RasGTP_parameter J_BRaf_iRaf1_bind J_neg_fdbk_phos_of_SOS J_neg_fdbk_unbind_piBRaf_RASt J_iBRaf_dim J_neg_fdbck_dephos_piRaf J_pRAF1_RAF1_phos J_neg_fdbk_unbind_RASt_nfpRAF1 ];
end

% -------------------------------------------------------
% ode rate
function dydt = f(t,y,p,y0)
	% State Variables
	Ras_iBRaf = y(1);
	Ras_Raf1_iBRaf_tetramer = y(2);
	Ras_BRaf = y(3);
	Ras_nfpRaf1 = y(4);
	iRaf1 = y(5);
	ERK = y(6);
	nfpSOS = y(7);
	Ras_pRaf1 = y(8);
	BRaf_iBRaf_dimer = y(9);
	Ras_pRaf1_iBRaf_tetramer = y(10);
	Ras_pRaf1_tetramer = y(11);
	Ras_nfpiRaf1 = y(12);
	Ras_GTP = y(13);
	nfpiRaf1 = y(14);
	Raf1 = y(15);
	Ras_BRaf_pRaf1_tetramer = y(16);
	Ras_pRaf1_Raf1_tetramer = y(17);
	nfpRaf1 = y(18);
	Ras_nfpBRaf = y(19);
	iBRaf = y(20);
	pMEK = y(21);
	mE = y(22);
	mEL = y(23);
	pRaf1 = y(24);
	EG2 = y(25);
	Ras_iRaf1 = y(26);
	Ras_nfpiBRaf = y(27);
	MEK = y(28);
	Ras_Braf_iRaf1_tetramer = y(29);
	Ras_Raf1 = y(30);
	mELmEL = y(31);
	nfpiBRaf = y(32);
	BRaf = y(33);
	SOS = y(34);
	nfpBRaf = y(35);
	GRB2_SOS = y(36);
	GRB2 = y(37);
	Ras_iRaf1_tetramer = y(38);
	iBRaf_dimer = y(39);
	Ras_GDP = y(40);
	E = y(41);
	Ras_BRaf_Raf1_tetramer = y(42);
	BRaf_dimer = y(43);
	EG2SOS = y(44);
	pERK = y(45);
	Ras_iRaf1_iBRaf1_tetramer = y(46);
	% Constants
	Ras_iRaf1_tetramer_init_molecules_um_2 = p(1);
	kEf = p(2);
	mlabfix_F_ = p(3);
	Ras_pRaf1_Raf1_tetramer_init_molecules_um_2 = p(4);
	kcatE = p(5);
	pRaf1_init_molecules_um_3 = p(6);
	netValence_neg_fdbk_phos_of_SOS = p(7);
	kpMEK = p(8);
	netValence_E_G2SOS_bind = p(9);
	knfpSOS = p(10);
	init_SOS_HeLa = p(11);
	netValence_pRaf1_Raf1_bind = p(12);
	netValence_iBRaf_dim = p(13);
	netValence_BRaf_iRaf1_bind = p(14);
	netValence_neg_fdbk_unbind_RASt_nfpRAF1 = p(15);
	kBr = p(16);
	Kr_ERK_phos = p(17);
	Kr_RAS_GTP_bind = p(18);
	init_ERK_HeLa = p(19);
	Kr_Ras_iRaf_nfp = p(20);
	netValence_pRAF1_RAF1_phos = p(21);
	kBf = p(22);
	kdpERK = p(23);
	Kr_pRAF1_dephos = p(24);
	kR1r = p(25);
	iRaf1_init_molecules_um_3 = p(26);
	kR1f = p(27);
	netValence_RAS_pRAF1_unbind = p(28);
	Kr_neg_fdbk_unbind_RASt_nfpBRAF = p(29);
	Kr_MEK_phos = p(30);
	init_MEK_HeLa = p(31);
	netValence_Ras_iBRaf_nfp = p(32);
	Kr_neg_fdbk_dephos_pRaf1 = p(33);
	netValence_RAS_iBRAF_bind = p(34);
	nfpiBRaf_init_molecules_um_3 = p(35);
	netValence_E_G2_bind = p(36);
	Kr_RAS_GTPhydro = p(37);
	Kr_Ras_iBRaf_nfp = p(38);
	pMEK_init_molecules_um_3 = p(39);
	Ras_BRaf_init_molecules_um_2 = p(40);
	Ras_iRaf1_iBRaf1_tetramer_init_molecules_um_2 = p(41);
	Ras_pRaf1_tetramer_init_molecules_um_2 = p(42);
	K_millivolts_per_volt = p(43);
	mlabfix_K_GHK_ = p(44);
	knfpBR1 = p(45);
	EG2_init_molecules_um_2 = p(46);
	Ras_BRaf_pRaf1_tetramer_init_molecules_um_2 = p(47);
	netValence_pRaf1_pRaf1_unbind = p(48);
	kpR1 = p(49);
	kSOSr = p(50);
	Ras_nfpRaf1_init_molecules_um_2 = p(51);
	netValence_iRaf1_iRaf1_bind = p(52);
	kSOSf = p(53);
	Ras_iBRaf_init_molecules_um_2 = p(54);
	netValence_RAS_BRAF_bind = p(55);
	Kr_ERK_dephos = p(56);
	kRgneslow = p(57);
	mlabfix_F_nmol_ = p(58);
	netValence_mEL_dim = p(59);
	netValence_BRAF_RAF1_bind = p(60);
	GRB2_SOS_init_molecules_um_3 = p(61);
	Kr_BRAF_RAF1_phos = p(62);
	E_init_molecules_um_2 = p(63);
	Ras_GTP_init_molecules_um_2 = p(64);
	kSon = p(65);
	kG2SOSr = p(66);
	kG2SOSf = p(67);
	Kr_neg_fdbk_dephos_iBRAF = p(68);
	nfpRaf1_init_molecules_um_3 = p(69);
	netValence_Raf1_iBRaf_phos = p(70);
	kiR1r = p(71);
	netValence_EG2_SOS_bind = p(72);
	netValence_neg_fdbk_unbind_RASt_nfpBRAF = p(73);
	mlabfix_T_ = p(74);
	kiR1f = p(75);
	Size_cytoplasm = p(76);
	kSoff = p(77);
	Kr_MEK_dephos = p(78);
	iBRaf_init_molecules_um_3 = p(79);
	netValence_BRaf_pRaf1_unbind = p(80);
	mlabfix_R_ = p(81);
	iBRaf_dimer_init_molecules_um_2 = p(82);
	netValence_BRAF_nfp = p(83);
	mELmEL_init_molecules_um_2 = p(84);
	Ras_Braf_iRaf1_tetramer_init_molecules_um_2 = p(85);
	kpERK = p(86);
	netValence_Ras_iRaf1_bind = p(87);
	Kr_neg_fdbk_unbind_piRaf1_RASt = p(88);
	netValence_Ras_pRaf1_dp = p(89);
	Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2 = p(90);
	netValence_EGF_bind = p(91);
	Kr_neg_fdbk_unbind_piBRaf_RASt = p(92);
	netValence_mELmEL_phos_dephos = p(93);
	EGF = p(94);
	kRhydro = p(95);
	Voltage_pm = p(96);
	mEL_init_molecules_um_2 = p(97);
	Kr_pRAF1_RAF1_phos = p(98);
	Kr_BRAF_nfp = p(99);
	netValence_pRaf1_iBRaf_unbind = p(100);
	netValence_Ras_Raf1_nfp = p(101);
	Kmgneslow = p(102);
	kdpMEK = p(103);
	knfpiR1r = p(104);
	kdR1r = p(105);
	kiBr = p(106);
	kdR1f = p(107);
	kiBf = p(108);
	init_GRB2_HeLa = p(109);
	Ras_nfpBRaf_init_molecules_um_2 = p(110);
	kfpBr = p(111);
	netValence_BRaf_iBRaf_dim = p(112);
	EG2SOS_init_molecules_um_2 = p(113);
	kdpSOS = p(114);
	init_RAF_HeLa = p(115);
	netValence_iRaf_iBRaf_dim = p(116);
	Kr_neg_fdbk_dephos_pSOS = p(117);
	netValence_neg_fdbk_unbind_piRaf1_RASt = p(118);
	netValence_iBRaf_Raf1_bind = p(119);
	kdp = p(120);
	netValence_neg_fdbk_unbind_piBRaf_RASt = p(121);
	netValence_RAS_GTP_bind = p(122);
	netValence_Ras_iRaf_nfp = p(123);
	nfpBRaf_init_molecules_um_3 = p(124);
	netValence_Ras_Raf1_bind = p(125);
	netValence_BRAF_RAF1_phos = p(126);
	kdEr = p(127);
	mlabfix_PI_ = p(128);
	kdEf = p(129);
	kdpR1 = p(130);
	init_BRAF_HeLa = p(131);
	BRaf_dimer_init_molecules_um_2 = p(132);
	Kr_neg_fdbck_dephos_piRaf = p(133);
	Kr_neg_fdbk_phos_of_SOS = p(134);
	Kr_neg_fdbk_unbind_RASt_nfpRAF1 = p(135);
	Ras_Raf1_iBRaf_tetramer_init_molecules_um_2 = p(136);
	BRaf_iBRaf_dimer_init_molecules_um_2 = p(137);
	init_RAS_HeLa = p(138);
	mlabfix_N_pmol_ = p(139);
	Ras_pRaf1_init_molecules_um_2 = p(140);
	init_EGFR_HeLa = p(141);
	nfpiRaf1_init_molecules_um_3 = p(142);
	kG2r = p(143);
	netValence_RAS_GTPhydro = p(144);
	Ras_nfpiRaf1_init_molecules_um_2 = p(145);
	init_sorafenib = p(146);
	kG2f = p(147);
	Ras_Raf1_init_molecules_um_2 = p(148);
	Kr_Ras_pRaf1_dp = p(149);
	KMOLE = p(150);
	Kr_Raf1_iBRaf_phos = p(151);
	Kr_neg_fdbk_dephos_pBRAF = p(152);
	pERK_init_molecules_um_3 = p(153);
	Ras_nfpiBRaf_init_molecules_um_2 = p(154);
	Size_pm = p(155);
	netValence_BRAF_dim = p(156);
	nfpSOS_init_molecules_um_3 = p(157);
	Ras_BRaf_Raf1_tetramer_init_molecules_um_2 = p(158);
	knfpiBr = p(159);
	Kr_Ras_Raf1_nfp = p(160);
	Ras_iRaf1_init_molecules_um_2 = p(161);
	kfpR1r = p(162);
	kEr = p(163);
	% Functions
	Kr_E_G2_bind = kG2r;
	Ras_GDP_init_molecules_um_2 = init_RAS_HeLa;
	Kf_iRaf_iBRaf_dim = kdR1f;
	Raf1_cyt_param = (Raf1 + pRaf1 + iRaf1 + nfpRaf1 + nfpiRaf1);
	Kr_E_G2SOS_bind = kG2SOSr;
	Kf_E_G2SOS_bind = kG2SOSf;
	J_E_G2SOS_bind = (((Kf_E_G2SOS_bind .* GRB2_SOS) .* E) - (Kr_E_G2SOS_bind .* EG2SOS));
	Kr_BRAF_dim = kdR1r;
	Kf_RAS_GTP_bind = (kRgneslow .* EG2SOS ./ (Kmgneslow + Ras_GDP));
	J_RAS_GTP_bind = ((Kf_RAS_GTP_bind .* Ras_GDP) - (Kr_RAS_GTP_bind .* Ras_GTP));
	Kf_Ras_iRaf_nfp = (kpR1 .* pERK);
	J_Ras_iRaf_nfp = ((Kf_Ras_iRaf_nfp .* Ras_iRaf1) - (Kr_Ras_iRaf_nfp .* Ras_nfpiRaf1));
	Kf_Ras_Raf1_nfp = (knfpBR1 .* pERK);
	Kf_Ras_iBRaf_nfp = (knfpBR1 .* pERK);
	Kf_pRAF1_dephos = kdpR1;
	J_pRAF1_dephos = ((Kf_pRAF1_dephos .* pRaf1) - (Kr_pRAF1_dephos .* Raf1));
	Kf_iRaf1_iRaf1_bind = kdR1f;
	Kf_BRAF_RAF1_phos = kpR1;
	Kf_BRAF_dim = kdR1f;
	J_BRAF_dim = ((Kf_BRAF_dim .* (Ras_BRaf ^ 2.0)) - (Kr_BRAF_dim .* BRaf_dimer));
	sorafenib_init_molecules_um_3 = init_sorafenib;
	sorafenib = sorafenib_init_molecules_um_3;
	Kf_E_G2_bind = kG2f;
	J_E_G2_bind = (((Kf_E_G2_bind .* GRB2) .* E) - (Kr_E_G2_bind .* EG2));
	Kf_BRAF_nfp = (knfpBR1 .* pERK);
	Kf_neg_fdbk_dephos_pSOS = kdpSOS;
	J_neg_fdbk_dephos_pSOS = ((Kf_neg_fdbk_dephos_pSOS .* nfpSOS) - (Kr_neg_fdbk_dephos_pSOS .* SOS));
	Kr_Ras_Raf1_bind = kR1r;
	Kf_neg_fdbk_dephos_iBRAF = kdpR1;
	Kf_ERK_dephos = kdpERK;
	Kf_iBRaf_dim = kdR1f;
	pRaf1_total = (Ras_pRaf1 + Ras_pRaf1_iBRaf_tetramer + Ras_pRaf1_tetramer + Ras_BRaf_pRaf1_tetramer + pRaf1 + Ras_pRaf1_Raf1_tetramer);
	Kf_EGF_bind = (kEf .* EGF);
	Kf_EG2_SOS_bind = kSOSf;
	Kf_RAS_GTPhydro = kRhydro;
	J_RAS_GTPhydro = ((Kf_RAS_GTPhydro .* Ras_GTP) - (Kr_RAS_GTPhydro .* Ras_GDP));
	Kf_iBRaf_Raf1_bind = kdR1f;
	Kf_neg_fdbk_unbind_RASt_nfpBRAF = kfpBr;
	Kr_mELmEL_phos_dephos = kdp;
	Kf_RAS_BRAF_bind = kBf;
	Kr_iBRaf_Raf1_bind = kdR1r;
	Kf_pRaf1_Raf1_bind = kdR1f;
	Kf_neg_fdbk_phos_of_SOS = (knfpSOS .* pERK);
	Kr_RAS_pRAF1_unbind = kR1f;
	Kf_neg_fdbk_unbind_piRaf1_RASt = knfpiR1r;
	Kf_Ras_Raf1_bind = kR1f;
	J_Ras_Raf1_bind = (((Kf_Ras_Raf1_bind .* Ras_GTP) .* Raf1) - (Kr_Ras_Raf1_bind .* Ras_Raf1));
	Kf_neg_fdbk_unbind_piBRaf_RASt = knfpiBr;
	Kf_Raf1_iBRaf_phos = kpR1;
	J_Raf1_iBRaf_phos = ((Kf_Raf1_iBRaf_phos .* Ras_Raf1_iBRaf_tetramer) - (Kr_Raf1_iBRaf_phos .* Ras_pRaf1_iBRaf_tetramer));
	Kr_BRaf_iBRaf_dim = kdR1r;
	Kf_BRaf_iBRaf_dim = kdR1f;
	Kf_sorafenib_Raf1_bind_cyt = kSon;
	J_ERK_dephos = ((Kf_ERK_dephos .* pERK) - (Kr_ERK_dephos .* ERK));
	Kf_Ras_pRaf1_dp = kdpR1;
	Kf_GRB2_SOS_bind_cyt = kSOSf;
	Kr_GRB2_SOS_bind_cyt = kSOSr;
	J_GRB2_SOS_bind_cyt = (((Kf_GRB2_SOS_bind_cyt .* GRB2) .* SOS) - (Kr_GRB2_SOS_bind_cyt .* GRB2_SOS));
	Kf_mEL_dim = kdEf;
	Kf_MEK_dephos = kdpMEK;
	J_neg_fdbk_unbind_RASt_nfpBRAF = ((Kf_neg_fdbk_unbind_RASt_nfpBRAF .* Ras_nfpBRaf) - ((Kr_neg_fdbk_unbind_RASt_nfpBRAF .* Ras_GTP) .* nfpBRaf));
	Kr_iRaf_iBRaf_dim = kdR1r;
	Kr_RAS_iBRAF_bind = kiBr;
	Kf_RAS_iBRAF_bind = kiBf;
	J_RAS_iBRAF_bind = (((Kf_RAS_iBRAF_bind .* iBRaf) .* Ras_GTP) - (Kr_RAS_iBRAF_bind .* Ras_iBRaf));
	Kf_neg_fdbk_dephos_pBRAF = kdpR1;
	J_neg_fdbk_dephos_pBRAF = ((Kf_neg_fdbk_dephos_pBRAF .* nfpBRaf) - (Kr_neg_fdbk_dephos_pBRAF .* BRaf));
	Kf_Ras_iRaf1_bind = kiR1f;
	Kr_Ras_iRaf1_bind = kiR1r;
	J_Ras_iRaf1_bind = (((Kf_Ras_iRaf1_bind .* iRaf1) .* Ras_GTP) - (Kr_Ras_iRaf1_bind .* Ras_iRaf1));
	SOS_init_molecules_um_3 = init_SOS_HeLa;
	mE_init_molecules_um_2 = init_EGFR_HeLa;
	KFlux_pm_cytoplasm = (Size_pm ./ Size_cytoplasm);
	ERK_init_molecules_um_3 = init_ERK_HeLa;
	Kr_EG2_SOS_bind = kSOSr;
	J_EG2_SOS_bind = (((Kf_EG2_SOS_bind .* EG2) .* SOS) - (Kr_EG2_SOS_bind .* EG2SOS));
	J_MEK_dephos = ((Kf_MEK_dephos .* pMEK) - (Kr_MEK_dephos .* MEK));
	Kf_ERK_phos = (kpERK .* pMEK);
	J_ERK_phos = ((Kf_ERK_phos .* ERK) - (Kr_ERK_phos .* pERK));
	Kf_BRAF_RAF1_bind = kdR1f;
	Kr_BRAF_RAF1_bind = kdR1r;
	J_BRAF_RAF1_bind = (((Kf_BRAF_RAF1_bind .* Ras_BRaf) .* Ras_Raf1) - (Kr_BRAF_RAF1_bind .* Ras_BRaf_Raf1_tetramer));
	Kf_pRaf1_pRaf1_unbind = kdR1r;
	Kr_pRaf1_pRaf1_unbind = kdR1f;
	J_pRaf1_pRaf1_unbind = ((Kf_pRaf1_pRaf1_unbind .* Ras_pRaf1_tetramer) - (Kr_pRaf1_pRaf1_unbind .* (Ras_pRaf1 ^ 2.0)));
	active_Raf = ((2.0 .* BRaf_dimer) + Ras_BRaf_Raf1_tetramer + (2.0 .* Ras_BRaf_pRaf1_tetramer) + Ras_pRaf1_iBRaf_tetramer + Ras_pRaf1 + Ras_pRaf1_Raf1_tetramer + (2.0 .* Ras_pRaf1_tetramer) + pRaf1 + BRaf_iBRaf_dimer + Ras_Braf_iRaf1_tetramer);
	Kf_MEK_phos = (kpMEK .* active_Raf);
	J_MEK_phos = ((Kf_MEK_phos .* MEK) - (Kr_MEK_phos .* pMEK));
	Kf_BRaf_pRaf1_unbind = kdR1r;
	Kr_BRaf_pRaf1_unbind = kdR1f;
	J_BRaf_pRaf1_unbind = ((Kf_BRaf_pRaf1_unbind .* Ras_BRaf_pRaf1_tetramer) - ((Kr_BRaf_pRaf1_unbind .* Ras_BRaf) .* Ras_pRaf1));
	Kf_BRaf_iRaf1_bind = kdR1f;
	Kf_pRaf1_iBRaf_unbind = kdR1r;
	Kr_pRaf1_iBRaf_unbind = kdR1f;
	J_pRaf1_iBRaf_unbind = ((Kf_pRaf1_iBRaf_unbind .* Ras_pRaf1_iBRaf_tetramer) - ((Kr_pRaf1_iBRaf_unbind .* Ras_iBRaf) .* Ras_pRaf1));
	Kf_sorafenib_BRaf_bind_cyt = kSon;
	J_iBRaf_Raf1_bind = (((Kf_iBRaf_Raf1_bind .* Ras_Raf1) .* Ras_iBRaf) - (Kr_iBRaf_Raf1_bind .* Ras_Raf1_iBRaf_tetramer));
	Kr_mEL_dim = kdEr;
	J_mEL_dim = ((Kf_mEL_dim .* (mEL ^ 2.0)) - (Kr_mEL_dim .* mELmEL));
	Kr_pRaf1_Raf1_bind = kdR1r;
	Kf_pRAF1_RAF1_phos = kpR1;
	Kf_RAS_pRAF1_unbind = kR1r;
	J_RAS_pRAF1_unbind = ((Kf_RAS_pRAF1_unbind .* Ras_pRaf1) - ((Kr_RAS_pRAF1_unbind .* pRaf1) .* Ras_GTP));
	Kr_BRaf_iRaf1_bind = kdR1r;
	Kf_neg_fdbck_dephos_piRaf = kdpR1;
	J_Ras_Raf1_nfp = ((Kf_Ras_Raf1_nfp .* Ras_Raf1) - (Kr_Ras_Raf1_nfp .* Ras_nfpRaf1));
	Kf_neg_fdbk_dephos_pRaf1 = kdpR1;
	J_neg_fdbk_dephos_pRaf1 = ((Kf_neg_fdbk_dephos_pRaf1 .* nfpRaf1) - (Kr_neg_fdbk_dephos_pRaf1 .* Raf1));
	Kr_sorafenib_Raf1_bind_cyt = kSoff;
	J_sorafenib_Raf1_bind_cyt = (((Kf_sorafenib_Raf1_bind_cyt .* sorafenib) .* Raf1) - (Kr_sorafenib_Raf1_bind_cyt .* iRaf1));
	Raf1_pm_param = (Ras_pRaf1 + (2.0 .* Ras_pRaf1_Raf1_tetramer) + Ras_BRaf_Raf1_tetramer + Ras_iRaf1 + Ras_Raf1 + (2.0 .* Ras_iRaf1_tetramer) + Ras_Braf_iRaf1_tetramer + Ras_BRaf_pRaf1_tetramer + (2.0 .* Ras_pRaf1_tetramer) + Ras_nfpiRaf1 + Ras_nfpRaf1 + Ras_iRaf1_iBRaf1_tetramer + Ras_Raf1_iBRaf_tetramer + Ras_pRaf1_iBRaf_tetramer);
	boundRaf_frac_parameter = (Raf1_pm_param ./ (Raf1_cyt_param + Raf1_pm_param));
	Kr_EGF_bind = kEr;
	J_BRAF_nfp = ((Kf_BRAF_nfp .* Ras_BRaf) - (Kr_BRAF_nfp .* Ras_nfpBRaf));
	J_Ras_iBRaf_nfp = ((Kf_Ras_iBRaf_nfp .* Ras_iBRaf) - (Kr_Ras_iBRaf_nfp .* Ras_nfpiBRaf));
	Kr_RAS_BRAF_bind = kBr;
	J_EGF_bind = ((Kf_EGF_bind .* mE) - (Kr_EGF_bind .* mEL));
	Kr_iRaf1_iRaf1_bind = kdR1r;
	MEK_init_molecules_um_3 = init_MEK_HeLa;
	J_iRaf1_iRaf1_bind = ((Kf_iRaf1_iRaf1_bind .* (Ras_iRaf1 ^ 2.0)) - (Kr_iRaf1_iRaf1_bind .* Ras_iRaf1_tetramer));
	GRB2_init_molecules_um_3 = init_GRB2_HeLa;
	J_BRaf_iBRaf_dim = (((Kf_BRaf_iBRaf_dim .* Ras_iBRaf) .* Ras_BRaf) - (Kr_BRaf_iBRaf_dim .* BRaf_iBRaf_dimer));
	Kf_mELmEL_phos_dephos = kcatE;
	J_mELmEL_phos_dephos = ((Kf_mELmEL_phos_dephos .* mELmEL) - (Kr_mELmEL_phos_dephos .* E));
	J_RAS_BRAF_bind = (((Kf_RAS_BRAF_bind .* BRaf) .* Ras_GTP) - (Kr_RAS_BRAF_bind .* Ras_BRaf));
	J_iRaf_iBRaf_dim = (((Kf_iRaf_iBRaf_dim .* Ras_iBRaf) .* Ras_iRaf1) - (Kr_iRaf_iBRaf_dim .* Ras_iRaf1_iBRaf1_tetramer));
	J_neg_fdbk_dephos_iBRAF = ((Kf_neg_fdbk_dephos_iBRAF .* nfpiBRaf) - (Kr_neg_fdbk_dephos_iBRAF .* iBRaf));
	J_Ras_pRaf1_dp = ((Kf_Ras_pRaf1_dp .* Ras_pRaf1) - (Kr_Ras_pRaf1_dp .* Ras_Raf1));
	Kr_iBRaf_dim = kdR1r;
	Kf_neg_fdbk_unbind_RASt_nfpRAF1 = kfpR1r;
	Raf1_init_molecules_um_3 = init_RAF_HeLa;
	J_BRAF_RAF1_phos = ((Kf_BRAF_RAF1_phos .* Ras_BRaf_Raf1_tetramer) - (Kr_BRAF_RAF1_phos .* Ras_BRaf_pRaf1_tetramer));
	Kr_sorafenib_BRaf_bind_cyt = kSoff;
	J_sorafenib_BRaf_bind_cyt = (((Kf_sorafenib_BRaf_bind_cyt .* sorafenib) .* BRaf) - (Kr_sorafenib_BRaf_bind_cyt .* iBRaf));
	J_pRaf1_Raf1_bind = (((Kf_pRaf1_Raf1_bind .* Ras_pRaf1) .* Ras_Raf1) - (Kr_pRaf1_Raf1_bind .* Ras_pRaf1_Raf1_tetramer));
	J_neg_fdbk_unbind_piRaf1_RASt = ((Kf_neg_fdbk_unbind_piRaf1_RASt .* Ras_nfpiRaf1) - ((Kr_neg_fdbk_unbind_piRaf1_RASt .* nfpiRaf1) .* Ras_GTP));
	BRaf_init_molecules_um_3 = init_BRAF_HeLa;
	RasGTP_parameter = (Ras_GTP + (2.0 .* (Ras_iRaf1_tetramer + Ras_Braf_iRaf1_tetramer + Ras_BRaf_Raf1_tetramer + Ras_pRaf1_Raf1_tetramer + Ras_BRaf_pRaf1_tetramer + Ras_pRaf1_tetramer + BRaf_dimer + BRaf_iBRaf_dimer + iBRaf_dimer + Ras_iRaf1_iBRaf1_tetramer + Ras_Raf1_iBRaf_tetramer + Ras_pRaf1_iBRaf_tetramer)) + Ras_BRaf + Ras_iRaf1 + Ras_pRaf1 + Ras_Raf1 + Ras_nfpRaf1 + Ras_nfpiRaf1 + Ras_iBRaf + Ras_nfpiBRaf + Ras_nfpBRaf);
	J_BRaf_iRaf1_bind = (((Kf_BRaf_iRaf1_bind .* Ras_BRaf) .* Ras_iRaf1) - (Kr_BRaf_iRaf1_bind .* Ras_Braf_iRaf1_tetramer));
	J_neg_fdbk_phos_of_SOS = ((Kf_neg_fdbk_phos_of_SOS .* EG2SOS) - ((Kr_neg_fdbk_phos_of_SOS .* EG2) .* nfpSOS));
	J_neg_fdbk_unbind_piBRaf_RASt = ((Kf_neg_fdbk_unbind_piBRaf_RASt .* Ras_nfpiBRaf) - ((Kr_neg_fdbk_unbind_piBRaf_RASt .* nfpiBRaf) .* Ras_GTP));
	J_iBRaf_dim = ((Kf_iBRaf_dim .* (Ras_iBRaf ^ 2.0)) - (Kr_iBRaf_dim .* iBRaf_dimer));
	J_neg_fdbck_dephos_piRaf = ((Kf_neg_fdbck_dephos_piRaf .* nfpiRaf1) - (Kr_neg_fdbck_dephos_piRaf .* iRaf1));
	J_pRAF1_RAF1_phos = ((Kf_pRAF1_RAF1_phos .* Ras_pRaf1_Raf1_tetramer) - (Kr_pRAF1_RAF1_phos .* Ras_pRaf1_tetramer));
	J_neg_fdbk_unbind_RASt_nfpRAF1 = ((Kf_neg_fdbk_unbind_RASt_nfpRAF1 .* Ras_nfpRaf1) - ((Kr_neg_fdbk_unbind_RASt_nfpRAF1 .* Ras_GTP) .* nfpRaf1));
	% Rates
	dydt = [
		(J_RAS_iBRAF_bind - (2.0 .* J_iBRaf_dim) - J_iBRaf_Raf1_bind + J_pRaf1_iBRaf_unbind - J_BRaf_iBRaf_dim - J_Ras_iBRaf_nfp - J_iRaf_iBRaf_dim);    % rate for Ras_iBRaf
		(J_iBRaf_Raf1_bind - J_Raf1_iBRaf_phos);    % rate for Ras_Raf1_iBRaf_tetramer
		( - J_BRAF_RAF1_bind + J_RAS_BRAF_bind - J_BRaf_iRaf1_bind + J_BRaf_pRaf1_unbind - (2.0 .* J_BRAF_dim) - J_BRAF_nfp - J_BRaf_iBRaf_dim);    % rate for Ras_BRaf
		(J_Ras_Raf1_nfp - J_neg_fdbk_unbind_RASt_nfpRAF1);    % rate for Ras_nfpRaf1
		(J_sorafenib_Raf1_bind_cyt - (KFlux_pm_cytoplasm .* J_Ras_iRaf1_bind) + J_neg_fdbck_dephos_piRaf);    % rate for iRaf1
		( - J_ERK_phos + J_ERK_dephos);    % rate for ERK
		( - J_neg_fdbk_dephos_pSOS + (KFlux_pm_cytoplasm .* J_neg_fdbk_phos_of_SOS));    % rate for nfpSOS
		( - J_RAS_pRAF1_unbind - J_pRaf1_Raf1_bind + (2.0 .* J_pRaf1_pRaf1_unbind) + J_BRaf_pRaf1_unbind + J_pRaf1_iBRaf_unbind - J_Ras_pRaf1_dp);    % rate for Ras_pRaf1
		J_BRaf_iBRaf_dim;    % rate for BRaf_iBRaf_dimer
		(J_Raf1_iBRaf_phos - J_pRaf1_iBRaf_unbind);    % rate for Ras_pRaf1_iBRaf_tetramer
		(J_pRAF1_RAF1_phos - J_pRaf1_pRaf1_unbind);    % rate for Ras_pRaf1_tetramer
		(J_Ras_iRaf_nfp - J_neg_fdbk_unbind_piRaf1_RASt);    % rate for Ras_nfpiRaf1
		( - J_Ras_Raf1_bind + J_RAS_pRAF1_unbind - J_RAS_BRAF_bind - J_Ras_iRaf1_bind - J_RAS_iBRAF_bind + J_neg_fdbk_unbind_RASt_nfpRAF1 + J_neg_fdbk_unbind_RASt_nfpBRAF + J_neg_fdbk_unbind_piRaf1_RASt + J_neg_fdbk_unbind_piBRaf_RASt + J_RAS_GTP_bind - J_RAS_GTPhydro);    % rate for Ras_GTP
		((KFlux_pm_cytoplasm .* J_neg_fdbk_unbind_piRaf1_RASt) - J_neg_fdbck_dephos_piRaf);    % rate for nfpiRaf1
		( - (KFlux_pm_cytoplasm .* J_Ras_Raf1_bind) - J_sorafenib_Raf1_bind_cyt + J_pRAF1_dephos + J_neg_fdbk_dephos_pRaf1);    % rate for Raf1
		(J_BRAF_RAF1_phos - J_BRaf_pRaf1_unbind);    % rate for Ras_BRaf_pRaf1_tetramer
		(J_pRaf1_Raf1_bind - J_pRAF1_RAF1_phos);    % rate for Ras_pRaf1_Raf1_tetramer
		((KFlux_pm_cytoplasm .* J_neg_fdbk_unbind_RASt_nfpRAF1) - J_neg_fdbk_dephos_pRaf1);    % rate for nfpRaf1
		(J_BRAF_nfp - J_neg_fdbk_unbind_RASt_nfpBRAF);    % rate for Ras_nfpBRaf
		(J_sorafenib_BRaf_bind_cyt - (KFlux_pm_cytoplasm .* J_RAS_iBRAF_bind) + J_neg_fdbk_dephos_iBRAF);    % rate for iBRaf
		(J_MEK_phos - J_MEK_dephos);    % rate for pMEK
		 - J_EGF_bind;    % rate for mE
		(J_EGF_bind - (2.0 .* J_mEL_dim));    % rate for mEL
		((KFlux_pm_cytoplasm .* J_RAS_pRAF1_unbind) - J_pRAF1_dephos);    % rate for pRaf1
		(J_E_G2_bind - J_EG2_SOS_bind + J_neg_fdbk_phos_of_SOS);    % rate for EG2
		(J_Ras_iRaf1_bind - (2.0 .* J_iRaf1_iRaf1_bind) - J_BRaf_iRaf1_bind - J_Ras_iRaf_nfp - J_iRaf_iBRaf_dim);    % rate for Ras_iRaf1
		(J_Ras_iBRaf_nfp - J_neg_fdbk_unbind_piBRaf_RASt);    % rate for Ras_nfpiBRaf
		( - J_MEK_phos + J_MEK_dephos);    % rate for MEK
		J_BRaf_iRaf1_bind;    % rate for Ras_Braf_iRaf1_tetramer
		(J_Ras_Raf1_bind - J_BRAF_RAF1_bind - J_pRaf1_Raf1_bind - J_iBRaf_Raf1_bind + J_Ras_pRaf1_dp - J_Ras_Raf1_nfp);    % rate for Ras_Raf1
		(J_mEL_dim - J_mELmEL_phos_dephos);    % rate for mELmEL
		((KFlux_pm_cytoplasm .* J_neg_fdbk_unbind_piBRaf_RASt) - J_neg_fdbk_dephos_iBRAF);    % rate for nfpiBRaf
		( - (KFlux_pm_cytoplasm .* J_RAS_BRAF_bind) - J_sorafenib_BRaf_bind_cyt + J_neg_fdbk_dephos_pBRAF);    % rate for BRaf
		( - J_GRB2_SOS_bind_cyt - (KFlux_pm_cytoplasm .* J_EG2_SOS_bind) + J_neg_fdbk_dephos_pSOS);    % rate for SOS
		((KFlux_pm_cytoplasm .* J_neg_fdbk_unbind_RASt_nfpBRAF) - J_neg_fdbk_dephos_pBRAF);    % rate for nfpBRaf
		( - (KFlux_pm_cytoplasm .* J_E_G2SOS_bind) + J_GRB2_SOS_bind_cyt);    % rate for GRB2_SOS
		( - (KFlux_pm_cytoplasm .* J_E_G2_bind) - J_GRB2_SOS_bind_cyt);    % rate for GRB2
		J_iRaf1_iRaf1_bind;    % rate for Ras_iRaf1_tetramer
		J_iBRaf_dim;    % rate for iBRaf_dimer
		( - J_RAS_GTP_bind + J_RAS_GTPhydro);    % rate for Ras_GDP
		(J_mELmEL_phos_dephos - J_E_G2SOS_bind - J_E_G2_bind);    % rate for E
		(J_BRAF_RAF1_bind - J_BRAF_RAF1_phos);    % rate for Ras_BRaf_Raf1_tetramer
		J_BRAF_dim;    % rate for BRaf_dimer
		(J_E_G2SOS_bind + J_EG2_SOS_bind - J_neg_fdbk_phos_of_SOS);    % rate for EG2SOS
		(J_ERK_phos - J_ERK_dephos);    % rate for pERK
		J_iRaf_iBRaf_dim;    % rate for Ras_iRaf1_iBRaf1_tetramer
	];
end
