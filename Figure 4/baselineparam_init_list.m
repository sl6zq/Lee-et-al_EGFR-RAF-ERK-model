
%from fullEGFR9_onemodel
params = [
	0.0;		% param(1) is 'Ras_iRaf1_tetramer_init_molecules_um_2'
	310.0;		% param(2) is 'kEf'
	96485.3321;		% param(3) is 'mlabfix_F_'
	0.0;		% param(4) is 'Ras_pRaf1_Raf1_tetramer_init_molecules_um_2'
	6.60624874;		% param(5) is 'kcatE'
	0.0;		% param(6) is 'pRaf1_init_molecules_um_3'
	1.0;		% param(7) is 'netValence_neg_fdbk_phos_of_SOS'
	6.0E-6;		% param(8) is 'kpMEK'
	1.0;		% param(9) is 'netValence_E_G2SOS_bind'
	0.0025;		% param(10) is 'knfpSOS'
	7700.0;		% param(11) is 'init_SOS_HeLa'
	1.0;		% param(12) is 'netValence_pRaf1_Raf1_bind'
	1.0;		% param(13) is 'netValence_iBRaf_dim'
	1.0;		% param(14) is 'netValence_BRaf_iRaf1_bind'
	1.0;		% param(15) is 'netValence_neg_fdbk_unbind_RASt_nfpRAF1'
	0.037;		% param(16) is 'kBr'
	0.0;		% param(17) is 'Kr_ERK_phos'
	0.0;		% param(18) is 'Kr_RAS_GTP_bind'
	600000.0;		% param(19) is 'init_ERK_HeLa'
	0.0;		% param(20) is 'Kr_Ras_iRaf_nfp'
	1.0;		% param(21) is 'netValence_pRAF1_RAF1_phos'
	2.5E-7;		% param(22) is 'kBf'
	0.6;		% param(23) is 'kdpERK'
	0.0;		% param(24) is 'Kr_pRAF1_dephos'
	8.4;		% param(25) is 'kR1r'
	0.0;		% param(26) is 'iRaf1_init_molecules_um_3'
	2.5E-5;		% param(27) is 'kR1f'
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
	6.0E-7;		% param(45) is 'knfpBR1'
	0.0;		% param(46) is 'EG2_init_molecules_um_2'
	0.0;		% param(47) is 'Ras_BRaf_pRaf1_tetramer_init_molecules_um_2'
	1.0;		% param(48) is 'netValence_pRaf1_pRaf1_unbind'
	190.0;		% param(49) is 'kpR1'
	0.0083;		% param(50) is 'kSOSr'
	0.0;		% param(51) is 'Ras_nfpRaf1_init_molecules_um_2'
	1.0;		% param(52) is 'netValence_iRaf1_iRaf1_bind'
	2.2E-6;		% param(53) is 'kSOSf'
	0.0;		% param(54) is 'Ras_iBRaf_init_molecules_um_2'
	1.0;		% param(55) is 'netValence_RAS_BRAF_bind'
	0.0;		% param(56) is 'Kr_ERK_dephos'
	5.838226569;		% param(57) is 'kRgneslow'
	9.64853321E-5;		% param(58) is 'mlabfix_F_nmol_'
	1.0;		% param(59) is 'netValence_mEL_dim'
	1.0;		% param(60) is 'netValence_BRAF_RAF1_bind'
	0.0;		% param(61) is 'GRB2_SOS_init_molecules_um_3'
	0.0;		% param(62) is 'Kr_BRAF_RAF1_phos'
	0.0;		% param(63) is 'E_init_molecules_um_2'
	0.0;		% param(64) is 'Ras_GTP_init_molecules_um_2'
	13.52143291;		% param(65) is 'kSon'
	4.6;		% param(66) is 'kG2SOSr'
	0.00942;		% param(67) is 'kG2SOSf'
	0.0;		% param(68) is 'Kr_neg_fdbk_dephos_iBRAF'
	0.0;		% param(69) is 'nfpRaf1_init_molecules_um_3'
	1.0;		% param(70) is 'netValence_Raf1_iBRaf_phos'
	0.072015388;		% param(71) is 'kiR1r'
	1.0;		% param(72) is 'netValence_EG2_SOS_bind'
	1.0;		% param(73) is 'netValence_neg_fdbk_unbind_RASt_nfpBRAF'
	300.0;		% param(74) is 'mlabfix_T_'
	1.15E-5;		% param(75) is 'kiR1f'
	1.0;		% param(76) is 'Size_cytoplasm'
	1243.069312;		% param(77) is 'kSoff'
	0.0;		% param(78) is 'Kr_MEK_dephos'
	0.0;		% param(79) is 'iBRaf_init_molecules_um_3'
	1.0;		% param(80) is 'netValence_BRaf_pRaf1_unbind'
	8314.46261815;		% param(81) is 'mlabfix_R_'
	0.0;		% param(82) is 'iBRaf_dimer_init_molecules_um_2'
	1.0;		% param(83) is 'netValence_BRAF_nfp'
	0.0;		% param(84) is 'mELmEL_init_molecules_um_2'
	0.0;		% param(85) is 'Ras_Braf_iRaf1_tetramer_init_molecules_um_2'
	6.0E-5;		% param(86) is 'kpERK'
	1.0;		% param(87) is 'netValence_Ras_iRaf1_bind'
	0.0;		% param(88) is 'Kr_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(89) is 'netValence_Ras_pRaf1_dp'
	0.0;		% param(90) is 'Ras_pRaf1_iBRaf_tetramer_init_molecules_um_2'
	1.0;		% param(91) is 'netValence_EGF_bind'
	0.0;		% param(92) is 'Kr_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(93) is 'netValence_mELmEL_phos_dephos'
	0.0017;		% param(94) is 'EGF'
	0.142907544;		% param(95) is 'kRhydro'
	0.0;		% param(96) is 'Voltage_pm'
	0.0;		% param(97) is 'mEL_init_molecules_um_2'
	0.0;		% param(98) is 'Kr_pRAF1_RAF1_phos'
	0.0;		% param(99) is 'Kr_BRAF_nfp'
	1.0;		% param(100) is 'netValence_pRaf1_iBRaf_unbind'
	1.0;		% param(101) is 'netValence_Ras_Raf1_nfp'
	141781.6814;		% param(102) is 'Kmgneslow'
	0.6;		% param(103) is 'kdpMEK'
	0.046229858;		% param(104) is 'knfpiR1r'
	110.0;		% param(105) is 'kdR1r'
	0.0367;		% param(106) is 'kiBr'
	0.019;		% param(107) is 'kdR1f'
	2.47E-7;		% param(108) is 'kiBf'
	628000.0;		% param(109) is 'init_GRB2_HeLa'
	0.0;		% param(110) is 'Ras_nfpBRaf_init_molecules_um_2'
	0.37;		% param(111) is 'kfpBr'
	1.0;		% param(112) is 'netValence_BRaf_iBRaf_dim'
	0.0;		% param(113) is 'EG2SOS_init_molecules_um_2'
	0.6;		% param(114) is 'kdpSOS'
	12000.0;		% param(115) is 'init_RAF_HeLa'
	1.0;		% param(116) is 'netValence_iRaf_iBRaf_dim'
	0.0;		% param(117) is 'Kr_neg_fdbk_dephos_pSOS'
	1.0;		% param(118) is 'netValence_neg_fdbk_unbind_piRaf1_RASt'
	1.0;		% param(119) is 'netValence_iBRaf_Raf1_bind'
	8.0;		% param(120) is 'kdp'
	1.0;		% param(121) is 'netValence_neg_fdbk_unbind_piBRaf_RASt'
	1.0;		% param(122) is 'netValence_RAS_GTP_bind'
	1.0;		% param(123) is 'netValence_Ras_iRaf_nfp'
	0.0;		% param(124) is 'nfpBRaf_init_molecules_um_3'
	1.0;		% param(125) is 'netValence_Ras_Raf1_bind'
	1.0;		% param(126) is 'netValence_BRAF_RAF1_phos'
	0.1;		% param(127) is 'kdEr'
	3.141592653589793;		% param(128) is 'mlabfix_PI_'
	9.2E-4;		% param(129) is 'kdEf'
	0.6;		% param(130) is 'kdpR1'
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
	460.0;		% param(143) is 'kG2r'
	1.0;		% param(144) is 'netValence_RAS_GTPhydro'
	0.0;		% param(145) is 'Ras_nfpiRaf1_init_molecules_um_2'
	10.0;		% param(146) is 'init_sorafenib'
	3.8E-4;		% param(147) is 'kG2f'
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
	0.37;		% param(159) is 'knfpiBr'
	0.0;		% param(160) is 'Kr_Ras_Raf1_nfp'
	0.0;		% param(161) is 'Ras_iRaf1_init_molecules_um_2'
	0.046229858;		% param(162) is 'kfpR1r'
	0.8;		% param(163) is 'kEr'
];



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

timeSpan=linspace(0,60);

%40 rate constants
wanted_param = [2; 5; 8; 10; 16; 22;
    23;
    25;
    27;
    45;
    49;
    50;
    53;
    57;
    65;
    66;
    67;
    71;
    75;
    77;
    86;
    95;
    102;
    103;
    104;
    105;
    106;
    107;
    108;
    111;
    114;
    120;
    127;
    129;
    130;
    143;
    147;
    159;
    162;
    163];


names = {'k_{E,f}'; 'k_{catE}'; 'k_{pMEK}'; 'k_{nfpSOS}'; 'k_{B,r}'; 'k_{B,f}'; 'k_{dpERK}'; 'k_{R1,r}'; 'k_{R1,f}'; 'k_{nfpBR1}'; 'k_{pR1}'; 'k_{SOS,r}'; 'k_{SOS,f}'; 'k_{Rgneslow}'; 'k_{Son}'; 'k_{G2SOS,r}'; 'k_{G2SOS,f}'; 'k_{iR1,r}'; 'k_{iR1,f}'; 'k_{Soff}'; 'k_{pERK}'; 'k_{Rhydro}';
'K_{mgneslow}'; 'k_{dpMEK}'; 'k_{nfpiR1r}'; 'k_{dR1,r}'; 'k_{iB,r}'; 'k_{dR1,f}'; 'k_{iB,f}'; 'k_{fpB,r}'; 'k_{dpSOS}';
'k_{dp}'; 'k_{dE,r}'; 'k_{dE,f}'; 'k_{dpR1}'; 'k_{G2,r}'; 'k_{G2,f}'; 'k_{nfpiB,r}'; 'k_{fpR1,r}'; 'k_{E,r}'};

PLSR_cat = categorical(names);
PLSR_cat = reordercats(PLSR_cat,{'k_{E,f}'; 'k_{catE}'; 'k_{pMEK}'; 'k_{nfpSOS}'; 'k_{B,r}'; 'k_{B,f}'; 'k_{dpERK}'; 'k_{R1,r}'; 'k_{R1,f}'; 'k_{nfpBR1}'; 'k_{pR1}'; 'k_{SOS,r}'; 'k_{SOS,f}'; 'k_{Rgneslow}'; 'k_{Son}'; 'k_{G2SOS,r}'; 'k_{G2SOS,f}'; 'k_{iR1,r}'; 'k_{iR1,f}'; 'k_{Soff}'; 'k_{pERK}'; 'k_{Rhydro}';
'K_{mgneslow}'; 'k_{dpMEK}'; 'k_{nfpiR1r}'; 'k_{dR1,r}'; 'k_{iB,r}'; 'k_{dR1,f}'; 'k_{iB,f}'; 'k_{fpB,r}'; 'k_{dpSOS}';
'k_{dp}'; 'k_{dE,r}'; 'k_{dE,f}'; 'k_{dpR1}'; 'k_{G2,r}'; 'k_{G2,f}'; 'k_{nfpiB,r}'; 'k_{fpR1,r}'; 'k_{E,r}'});
