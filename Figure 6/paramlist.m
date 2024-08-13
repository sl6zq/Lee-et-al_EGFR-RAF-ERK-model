% Best-fit parameter values from model fitting in SloppyCell - conducted by co-PI Kevin Brown 
% Define initial conditions for ODE model, indices and names of parameters of interest

params = [
	0.0;		% param(1) is 'Ras_iRaf1_tetramer_init_molecules_um_2'
	1.251011e+04;		% param(2) is 'kEf'
	96485.3321;		% param(3) is 'mlabfix_F_'
	0.0;		% param(4) is 'Ras_pRaf1_Raf1_tetramer_init_molecules_um_2'
	6.663433e+00;		% param(5) is 'kcatE'
	0.0;		% param(6) is 'pRaf1_init_molecules_um_3'
	1.0;		% param(7) is 'netValence_neg_fdbk_phos_of_SOS'
	1.587740e-04;	% param(8) is 'kpMEK'
	1.0;		% param(9) is 'netValence_E_G2SOS_bind'
	7.836388e-06;		% param(10) is 'knfpSOS'
	7700.0;		% param(11) is 'init_SOS_HeLa'
	1.0;		% param(12) is 'netValence_pRaf1_Raf1_bind'
	1.0;		% param(13) is 'netValence_iBRaf_dim'
	1.0;		% param(14) is 'netValence_BRaf_iRaf1_bind'
	1.0;		% param(15) is 'netValence_neg_fdbk_unbind_RASt_nfpRAF1'
	1.726000e-02;		% param(16) is 'kBr'
	0.0;		% param(17) is 'Kr_ERK_phos'
	0.0;		% param(18) is 'Kr_RAS_GTP_bind'
	600000.0;		% param(19) is 'init_ERK_HeLa'
	0.0;		% param(20) is 'Kr_Ras_iRaf_nfp'
	1.0;		% param(21) is 'netValence_pRAF1_RAF1_phos'
	5.634366e+02;		% param(22) is 'kBf'
	3.766460e-01;		% param(23) is 'kdpERK'
	0.0;		% param(24) is 'Kr_pRAF1_dephos'
	9.084007e+01;	% param(25) is 'kR1r'
	0.0;		% param(26) is 'iRaf1_init_molecules_um_3'
	3.314178e-05;		% param(27) is 'kR1f'
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
	3.805618e-13;		% param(45) is 'knfpBR1'
	0.0;		% param(46) is 'EG2_init_molecules_um_2'
	0.0;		% param(47) is 'Ras_BRaf_pRaf1_tetramer_init_molecules_um_2'
	1.0;		% param(48) is 'netValence_pRaf1_pRaf1_unbind'
	2.288588e+03;		% param(49) is 'kpR1'
	1.905146e+00;		% param(50) is 'kSOSr'
	0.0;		% param(51) is 'Ras_nfpRaf1_init_molecules_um_2'
	1.0;		% param(52) is 'netValence_iRaf1_iRaf1_bind'
	5.520185e-05;		% param(53) is 'kSOSf'
	0.0;		% param(54) is 'Ras_iBRaf_init_molecules_um_2'
	1.0;		% param(55) is 'netValence_RAS_BRAF_bind'
	0.0;		% param(56) is 'Kr_ERK_dephos'
	5.754841e+00;		% param(57) is 'kRgneslow'
	9.64853321E-5;		% param(58) is 'mlabfix_F_nmol_'
	1.0;		% param(59) is 'netValence_mEL_dim'
	1.0;		% param(60) is 'netValence_BRAF_RAF1_bind'
	0.0;		% param(61) is 'GRB2_SOS_init_molecules_um_3'
	0.0;		% param(62) is 'Kr_BRAF_RAF1_phos'
	0.0;		% param(63) is 'E_init_molecules_um_2'
	0.0;		% param(64) is 'Ras_GTP_init_molecules_um_2'
	2.828193e+01;		% param(65) is 'kSon'
	5.887828e+01;		% param(66) is 'kG2SOSr'
	2.594363e-01;    	% param(67) is 'kG2SOSf'
	0.0;		% param(68) is 'Kr_neg_fdbk_dephos_iBRAF'
	0.0;		% param(69) is 'nfpRaf1_init_molecules_um_3'
	1.0;		% param(70) is 'netValence_Raf1_iBRaf_phos'
	2.201878e-01;		% param(71) is 'kiR1r'
	1.0;		% param(72) is 'netValence_EG2_SOS_bind'
	1.0;		% param(73) is 'netValence_neg_fdbk_unbind_RASt_nfpBRAF'
	300.0;		% param(74) is 'mlabfix_T_'
	1.480046e-06;		% param(75) is 'kiR1f'
	1.0;		% param(76) is 'Size_cytoplasm'
	8.722582e+00;		% param(77) is 'kSoff'
	0.0;		% param(78) is 'Kr_MEK_dephos'
	0.0;		% param(79) is 'iBRaf_init_molecules_um_3'
	1.0;		% param(80) is 'netValence_BRaf_pRaf1_unbind'
	8314.46261815;		% param(81) is 'mlabfix_R_'
	0.0;		% param(82) is 'iBRaf_dimer_init_molecules_um_2'
	1.0;		% param(83) is 'netValence_BRAF_nfp'
	0.0;		% param(84) is 'mELmEL_init_molecules_um_2'
	0.0;		% param(85) is 'Ras_Braf_iRaf1_tetramer_init_molecules_um_2'
	1.161078e-06;		% param(86) is 'kpERK'
	1.0;		% param(87) is 'netValence_Ras_iRaf1_bind'
	0.0;		% param(88) is 'Kr_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(89) is 'netValence_Ras_pRaf1_dp'
	0.0;		% param(90) is 'Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2'
	1.0;		% param(91) is 'netValence_EGF_bind'
	0.0;		% param(92) is 'Kr_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(93) is 'netValence_mELmEL_phos_dephos'
	0.0017;		% param(94) is 'EGF'
	2.711843e-01;		% param(95) is 'kRhydro'
	0.0;		% param(96) is 'Voltage_pm'
	0.0;		% param(97) is 'mEL_init_molecules_um_2'
	0.0;		% param(98) is 'Kr_pRAF1_RAF1_phos'
	0.0;		% param(99) is 'Kr_BRAF_nfp'
	1.0;		% param(100) is 'netValence_pRaf1_iBRaf_unbind'
	1.0;		% param(101) is 'netValence_Ras_Raf1_nfp'
	7.205058e+04;		% param(102) is 'Kmgneslow'
	5.109929e-01;		% param(103) is 'kdpMEK'
	1.741575e-02;		% param(104) is 'knfpiR1r'
	1.314078e+01; % param(105) is 'kdR1r'
	1.564178e+00;		% param(106) is 'kiBr'
	9.753946e-01; 	% param(107) is 'kdR1f'
	1.579755e-13;		% param(108) is 'kiBf'
	628000.0;		% param(109) is 'init_GRB2_HeLa'
	0.0;		% param(110) is 'Ras_nfpBRaf_init_molecules_um_2'
	5.322688e+00; 	% param(111) is 'kfpBr'
	1.0;		% param(112) is 'netValence_BRaf_iBRaf_dim'
	0.0;		% param(113) is 'EG2SOS_init_molecules_um_2'
	5.861134e-03; 		% param(114) is 'kdpSOS'
	12000.0;		% param(115) is 'init_RAF_HeLa'
	1.0;		% param(116) is 'netValence_iRaf_iBRaf_dim'
	0.0;		% param(117) is 'Kr_neg_fdbk_dephos_pSOS'
	1.0;		% param(118) is 'netValence_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(119) is 'netValence_iBRaf_Raf1_bind'
	5.178373e+01;		% param(120) is 'kdp'
	1.0;		% param(121) is 'netValence_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(122) is 'netValence_RAS_GTP_bind'
	1.0;		% param(123) is 'netValence_Ras_iRaf_nfp'
	0.0;		% param(124) is 'nfpBRaf_init_molecules_um_3'
	1.0;		% param(125) is 'netValence_Ras_Raf1_bind'
	1.0;		% param(126) is 'netValence_BRAF_RAF1_phos'
	1.518401e+01;		% param(127) is 'kdEr'
	3.141592653589793;		% param(128) is 'mlabfix_PI_'
	9.793454e-01;		% param(129) is 'kdEf'
	2.812755e+01; 		% param(130) is 'kdpR1'
	1000.0;		% param(131) is 'init_BRAF_HeLa'
	0.0;		% param(132) is 'BRaf_dimer_init_molecules_um_2'
	0.0;		% param(133) is 'Kr_neg_fdbck_dephos_piRaf'
	0.0;		% param(134) is 'Kr_neg_fdbk_phos_of_SOS'
	0.0;		% param(135) is 'Kr_neg_fdbk_unbind_RASt_nfpRAF1'
	0.0;		% param(136) is 'Ras_Raf1_iBRaf_tetramer_init_molecules_um_2'
	0.0;		% param(137) is 'BRaf_iBRaf_dimer_init_molecules_um_2'
	135000.0;		% param(138) is 'init_RAS_HeLa'
	6.02214179E11;		% param(139) is 'mlabfix_N_pmol_'
	0.0;		% param(140) is 'Ras_pRaf1_init_molecules_um_2'
	93000.0;		% param(141) is 'init_EGFR_HeLa'
	0.0;		% param(142) is 'nfpiRaf1_init_molecules_um_3'
	3.333723e+03;	% param(143) is 'kG2r'
	1.0;		% param(144) is 'netValence_RAS_GTPhydro'
	0.0;		% param(145) is 'Ras_nfpiRaf1_init_molecules_um_2'
	10.0;		% param(146) is 'init_sorafenib'
	6.520470e-01;		% param(147) is 'kG2f'
	0.0;		% param(148) is 'Ras_Raf1_init_molecules_um_2'
	0.0;		% param(149) is 'Kr_Ras_pRaf1_dp'
	0.001660538783162726;		% param(150) is 'KMOLE'
	0.0;		% param(151) is 'Kr_Raf1_iBRaf_phos'
	0.0;		% param(152) is 'Kr_neg_fdbk_dephos_pBRAF'
	0.0;		% param(153) is 'pERK_init_molecules_um_3'
	0.0;		% param(154) is 'Ras_nfpiBRaf_init_molecules_um_2'
	1.0;		% param(155) is 'Size_pm'
	1.0;		% param(156) is 'netValence_BRAF_dim'
	0.0;		% param(157) is 'nfpSOS_init_molecules_um_3'
	0.0;		% param(158) is 'Ras_BRaf_Raf1_tetramer_init_molecules_um_2'
	1.022145e+01;		% param(159) is 'knfpiBr'
	0.0;		% param(160) is 'Kr_Ras_Raf1_nfp'
	0.0;		% param(161) is 'Ras_iRaf1_init_molecules_um_2'
	4.116000e-02;		% param(162) is 'kfpR1r'
	2.157247e+01; 		% param(163) is 'kEr'
];

% Initial conditions for ODE model

Ras_GDP_init_molecules_um_2 = params(138);
ERK_init_molecules_um_3 = params(19);
Raf1_init_molecules_um_3 = params(115);
mE_init_molecules_um_2 = params(141);
MEK_init_molecules_um_3 = params(31);
BRaf_init_molecules_um_3 = params(131);
SOS_init_molecules_um_3 = params(11);
GRB2_init_molecules_um_3 = params(109);	

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

% Indices and names of parameters of interest
% 40 rate constants (binding on & off)
wanted_param = [2; 5; 8; 10; 16; 22; 23; 25; 27; 45; ...
    49; 50; 53; 57; 65; 66; 67; 71; 75; 77; ...
    86; 95; 102; 103; 104; 105; 106; 107; 108; 111; ...
    114; 120; 127; 129; 130; 143; 147; 159; 162; 163];

names = {'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdE,f'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
names_vec = ["kEf" "kcatE" "kpMEK" "knfpSOS" "kBr" "kBf" "kdpERK" "kR1r" "kR1f" "knfpBR1" "kpR1" "kSOSr" "kSOSf" "kRgneslow" "kSon" "kG2SOSr" "kG2SOSf" "kiR1r" "kiR1f" "kSoff" "kpERK" "kRhydro" "Kmgneslow" "kdpMEK" "knfpiR1r" "kdR1r" "kiBr" "kdR1f" "kiBf" "kfpBr" "kdpSOS" "kdp" "kdEr" "kdE,f" "kdpR1" "kG2r" "kG2f" "knfpiBr" "kfpR1r" "kEr"]';
PLSR_cat = categorical(names);
PLSR_cat = reordercats(PLSR_cat,{'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdE,f'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'});

names_porder    = {'kEf'; 'kEr'; 'kdEf'; 'kdEr';...
    'kcatE'; 'kdp';...
    'kG2f'; 'kG2r'; 'kSOSf'; 'kSOSr';'kG2SOSf';'kG2SOSr'; 'knfpSOS'; 'kdpSOS';...
    'kRgneslow'; 'Kmgneslow';'kRhydro'; ...
    'kR1f'; 'kR1r'; 'kdR1f'; 'kdR1r'; 'kfpR1r'; 'kpR1'; 'kdpR1';...
    'kBf'; 'kBr'; 'knfpBR1'; 'kfpBr'; 'knfpiBr';...
    'kSon';  'kSoff';  'kiR1f'; 'kiR1r'; 'kiBf'; 'kiBr';  'knfpiR1r'; ...
    'kpMEK';'kdpMEK';'kpERK';'kdpERK'};
PLSR_cat_porder = categorical(names_porder);
PLSR_cat_porder = reordercats(PLSR_cat_porder,{'kEf'; 'kEr'; 'kdEf'; 'kdEr';...
    'kcatE'; 'kdp';...
    'kG2f'; 'kG2r'; 'kSOSf'; 'kSOSr';'kG2SOSf';'kG2SOSr'; 'knfpSOS'; 'kdpSOS';...
    'kRgneslow'; 'Kmgneslow';'kRhydro'; ...
    'kR1f'; 'kR1r'; 'kdR1f'; 'kdR1r'; 'kfpR1r'; 'kpR1'; 'kdpR1';...
    'kBf'; 'kBr'; 'knfpBR1'; 'kfpBr'; 'knfpiBr';...
    'kSon';  'kSoff';  'kiR1f'; 'kiR1r'; 'kiBf'; 'kiBr';  'knfpiR1r'; ...
    'kpMEK';'kdpMEK';'kpERK';'kdpERK'});
wanted_param_porder = [2;163;129;127;5;120;147;143;53;50;67;66;10;114;...
57;102;95;27;25;107;105;162;49;130;22;16;45;111;159;65;77;75;71;108;106;...
104;8;103;86;23];
% rate constants (off only)
wanted_param_offonly = [16; 23; 25; 50; 66; 71; 77; 86; 95; 103; ...
    104; 105; 106; 114; 120; 127; 143; 159; 162; 163];

names_offonly = {'kBr'; 'kdpERK'; 'kR1r'; 'kSOSr'; 'kG2SOSr'; 'kiR1r'; 'kSoff'; 'kRhydro';
'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kfpBr'; 'kdpSOS';'kdp';
'kdEr'; 'kG2r'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
PLSR_cat_offonly = categorical(names_offonly);
PLSR_cat_offonly = reordercats(PLSR_cat_offonly,{'kBr'; 'kdpERK'; 'kR1r'; 'kSOSr'; 'kG2SOSr'; 'kiR1r'; 'kSoff'; 'kRhydro';
'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kfpBr'; 'kdpSOS';'kdp';
'kdEr'; 'kG2r'; 'knfpiBr'; 'kfpR1r'; 'kEr'});

% Expression levels of key proteins + EGF concentration  
wanted_param_exp = [94; 141; 109; 11; 138; 115; 131; 31; 19]; %94: EGF, 141: EGFR, 109: GRB2, 11: SOS, 138: RAS, 115: RAF, 131: BRAF, 31: MEK, 19: ERK, 146: sorafenib
names_exp = {'[EGF]';'[EGFR]';'[GRB2]';'[SOS]';'[Ras-GDP]';'[Raf1]'; '[BRaf]'; '[MEK]'; '[ERK]'};


names_all = {'Ras_iRaf1_tetramer_init_molecules_um_2';
	'kEf';
	'mlabfix_F_';
	'Ras_pRaf1_Raf1_tetramer_init_molecules_um_2';
	'kcatE';
	'pRaf1_init_molecules_um_3';
	'netValence_neg_fdbk_phos_of_SOS';
	'kpMEK';
	'netValence_E_G2SOS_bind';
	'knfpSOS';
	'init_SOS_HeLa';
	'netValence_pRaf1_Raf1_bind';
	'netValence_iBRaf_dim';
	'netValence_BRaf_iRaf1_bind';
	'netValence_neg_fdbk_unbind_RASt_nfpRAF1';
	'kBr';
	'Kr_ERK_phos';
	'Kr_RAS_GTP_bind';
	'init_ERK_HeLa';
	'Kr_Ras_iRaf_nfp';
	'netValence_pRAF1_RAF1_phos';
	'kBf';
	'kdpERK';
	'Kr_pRAF1_dephos';
	'kR1r';
	'iRaf1_init_molecules_um_3';
	'kR1f';
	'netValence_RAS_pRAF1_unbind';
	'Kr_neg_fdbk_unbind_RASt_nfpBRAF';
	'Kr_MEK_phos';
	'init_MEK_HeLa';
	'netValence_Ras_iBRaf_nfp';
	'Kr_neg_fdbk_dephos_pRaf1';
	'netValence_RAS_iBRAF_bind';
	'nfpiBRaf_init_molecules_um_3';
	'netValence_E_G2_bind';
	'Kr_RAS_GTPhydro';
	'Kr_Ras_iBRaf_nfp';
	'pMEK_init_molecules_um_3';
	'Ras_BRaf_init_molecules_um_2';
	'Ras_iRaf1_iBRaf1_tetramer_init_molecules_um_2';
	'Ras_pRaf1_tetramer_init_molecules_um_2';
	'K_millivolts_per_volt';
	'mlabfix_K_GHK_';
	'knfpBR1';
	'EG2_init_molecules_um_2';
	'Ras_BRaf_pRaf1_tetramer_init_molecules_um_2';
	'netValence_pRaf1_pRaf1_unbind';
	'kpR1';
	'kSOSr';
	'Ras_nfpRaf1_init_molecules_um_2';
	'netValence_iRaf1_iRaf1_bind';
	'kSOSf';
	'Ras_iBRaf_init_molecules_um_2';
	'netValence_RAS_BRAF_bind';
	'Kr_ERK_dephos';
	'kRgneslow';
	'mlabfix_F_nmol_';
	'netValence_mEL_dim';
	'netValence_BRAF_RAF1_bind';
	'GRB2_SOS_init_molecules_um_3';
	'Kr_BRAF_RAF1_phos';
	'E_init_molecules_um_2';
	'Ras_GTP_init_molecules_um_2';
	'kSon';
	'kG2SOSr';
	'kG2SOSf';
	'Kr_neg_fdbk_dephos_iBRAF';
	'nfpRaf1_init_molecules_um_3';
	'netValence_Raf1_iBRaf_phos';
	'kiR1r';
	'netValence_EG2_SOS_bind';
	'netValence_neg_fdbk_unbind_RASt_nfpBRAF';
	'mlabfix_T_';
	'kiR1f';
	'Size_cytoplasm';
	'kSoff';
	'Kr_MEK_dephos';
	'iBRaf_init_molecules_um_3';
	'netValence_BRaf_pRaf1_unbind';
	'mlabfix_R_';
	'iBRaf_dimer_init_molecules_um_2';
	'netValence_BRAF_nfp';
	'mELmEL_init_molecules_um_2';
	'Ras_Braf_iRaf1_tetramer_init_molecules_um_2';
	'kpERK';
	'netValence_Ras_iRaf1_bind';
	'Kr_neg_fdbk_unbind_piRaf1_RASt';
	'netValence_Ras_pRaf1_dp';
	'Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2';
	'netValence_EGF_bind';
	'Kr_neg_fdbk_unbind_piBRaf_RASt';
	'netValence_mELmEL_phos_dephos';
	'EGF';
	'kRhydro';
	'Voltage_pm';
	'mEL_init_molecules_um_2';
	'Kr_pRAF1_RAF1_phos';
	'Kr_BRAF_nfp';
	'netValence_pRaf1_iBRaf_unbind';
	'netValence_Ras_Raf1_nfp';
	'Kmgneslow';
	'kdpMEK';
	'knfpiR1r';
	'kdR1r';
	'kiBr';
	'kdR1f';
	'kiBf';
	'init_GRB2_HeLa';
	'Ras_nfpBRaf_init_molecules_um_2';
	'kfpBr';
	'netValence_BRaf_iBRaf_dim';
	'EG2SOS_init_molecules_um_2';
	'kdpSOS';
	'init_RAF_HeLa';
	'netValence_iRaf_iBRaf_dim';
	'Kr_neg_fdbk_dephos_pSOS';
	'netValence_neg_fdbk_unbind_piRaf1_RASt';
	'netValence_iBRaf_Raf1_bind';
	'kdp';
	'netValence_neg_fdbk_unbind_piBRaf_RASt';
	'netValence_RAS_GTP_bind';
	'netValence_Ras_iRaf_nfp';
	'nfpBRaf_init_molecules_um_3';
	'netValence_Ras_Raf1_bind';
	'netValence_BRAF_RAF1_phos';
	'kdEr';
	'mlabfix_PI_';
	'kdEf';
	'kdpR1';
	'init_BRAF_HeLa';
	'BRaf_dimer_init_molecules_um_2';
	'Kr_neg_fdbck_dephos_piRaf';
	'Kr_neg_fdbk_phos_of_SOS';
	'Kr_neg_fdbk_unbind_RASt_nfpRAF1';
	'Ras_Raf1_iBRaf_tetramer_init_molecules_um_2';
	'BRaf_iBRaf_dimer_init_molecules_um_2';
	'init_RAS_HeLa';
	'mlabfix_N_pmol_';
	'Ras_pRaf1_init_molecules_um_2';
	'init_EGFR_HeLa';
	'nfpiRaf1_init_molecules_um_3';
	'kG2r';
	'netValence_RAS_GTPhydro';
	'Ras_nfpiRaf1_init_molecules_um_2';
	'init_sorafenib';
	'kG2f';
	'Ras_Raf1_init_molecules_um_2';
	 'Kr_Ras_pRaf1_dp';
	'KMOLE';
	'Kr_Raf1_iBRaf_phos';
	'Kr_neg_fdbk_dephos_pBRAF';
	'pERK_init_molecules_um_3';
	'Ras_nfpiBRaf_init_molecules_um_2';
	'Size_pm';
	'netValence_BRAF_dim';
	'nfpSOS_init_molecules_um_3';
	'Ras_BRaf_Raf1_tetramer_init_molecules_um_2';
	'knfpiBr';
	'Kr_Ras_Raf1_nfp';
	'Ras_iRaf1_init_molecules_um_2';
	'kfpR1r';
	'kEr';
};