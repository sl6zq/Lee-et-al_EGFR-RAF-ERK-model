%Select only rate constants (does not include kRhydro)
wanted_param_mutant                = [2; 5; 8; 10; 16; 22; 23; 25; 27; 45;49;50;53;57;65;66; 67;71;75;77;86;102;103;104;105;106;107;108;111;114;120;127;129;130;143;147;159;162;163];

wanted_param_wt                    = wanted_param_mutant;

names_mutant                       = {'\itk{E,f}'; '\itk{catE}'; '\itk{pMEK}'; '\itk{nfpSOS}'; ...
                                      '\itk{B,r}'; '\itk{B,f}'; '\itk{dpERK}'; '\itk{R1,r}'; '\itk{R1,f}'; '\itk{nfpBR1}'; ...
                                      '\itk{pR1}'; '\itk{SOS,r}'; '\itk{SOS,f}'; '\itk{Rgneslow}'; '\itk{Son}'; '\itk{G2SOS,r}'; ...
                                      '\itk{G2SOS,f}'; '\itk{iR1,r}'; '\itk{iR1,f}'; '\itk{Soff}'; '\itk{pERK}';...
                                      '\itK{mgneslow}'; '\itk{dpMEK}'; '\itk{nfpiR1r}'; '\itk{dR1,r}'; '\itk{iB,r}'; ...
                                      '\itk{dR1,f}'; '\itk{iB,f}'; '\itk{fpB,r}'; '\itk{dpSOS}';'\itk{dp}'; '\itk{dE,r}';...
                                      '\itk{dE,f}'; '\itk{dpR1}'; '\itk{G2,r}'; '\itk{G2,f}'; '\itk{nfpiB,r}'; '\itk{fpR1,r}'; '\itk{E,r}'};
PLSR_cat_mutant                    = categorical(names_mutant);
PLSR_cat_mutant                    = reordercats(PLSR_cat_mutant,{'\itk{E,f}'; '\itk{catE}'; '\itk{pMEK}'; '\itk{nfpSOS}';...
                                     '\itk{B,r}'; '\itk{B,f}'; '\itk{dpERK}'; '\itk{R1,r}'; '\itk{R1,f}'; '\itk{nfpBR1}'; '\itk{pR1}'; '\itk{SOS,r}'; ...
                                     '\itk{SOS,f}'; '\itk{Rgneslow}'; '\itk{Son}'; '\itk{G2SOS,r}'; '\itk{G2SOS,f}'; '\itk{iR1,r}'; '\itk{iR1,f}'; ...
                                     '\itk{Soff}'; '\itk{pERK}';'\itK{mgneslow}'; '\itk{dpMEK}'; '\itk{nfpiR1r}'; '\itk{dR1,r}'; '\itk{iB,r}'; ...
                                     '\itk{dR1,f}'; '\itk{iB,f}'; '\itk{fpB,r}'; '\itk{dpSOS}'; '\itk{dp}'; '\itk{dE,r}'; '\itk{dE,f}'; '\itk{dpR1}';...
                                     '\itk{G2,r}'; '\itk{G2,f}'; '\itk{nfpiB,r}'; '\itk{fpR1,r}'; '\itk{E,r}'});

%kRhydro is 16-fold lower in Ras mutants
params_mutant                      = params;
params_mutant(95)                  = params_mutant(95)/16; %KRASG12V kRhydro

%Define initial conditions
Ras_GDP_init_molecules_um_2_mutant = params_mutant(138);
ERK_init_molecules_um_3_mutant     = params_mutant(19);
Raf1_init_molecules_um_3_mutant    = params_mutant(115);
mE_init_molecules_um_2_mutant      = params_mutant(141);
MEK_init_molecules_um_3_mutant     = params_mutant(31);
BRaf_init_molecules_um_3_mutant    = params_mutant(131);
SOS_init_molecules_um_3_mutant     = params_mutant(11);
GRB2_init_molecules_um_3_mutant    = params_mutant(109);	

yinit_mutant                       = [
                                      0.0;		% yinit(1) is the initial condition for 'Ras_iBRaf'
                                      0.0;		% yinit(2) is the initial condition for 'Ras_Raf1_iBRaf_tetramer'
                                      0.0;		% yinit(3) is the initial condition for 'Ras_BRaf'
                                      0.0;		% yinit(4) is the initial condition for 'Ras_nfpRaf1'
                                      0.0;		% yinit(5) is the initial condition for 'iRaf1'
                                      ERK_init_molecules_um_3_mutant;		% yinit(6) is the initial condition for 'ERK'
                                      0.0;		% yinit(7) is the initial condition for 'nfpSOS'
                                      0.0;		% yinit(8) is the initial condition for 'Ras_pRaf1'
                                      0.0;		% yinit(9) is the initial condition for 'BRaf_iBRaf_dimer'
                                      0.0;		% yinit(10) is the initial condition for 'Ras_pRaf1_iBRaf_tetramer'
                                      0.0;		% yinit(11) is the initial condition for 'Ras_pRaf1_tetramer'
                                      0.0;		% yinit(12) is the initial condition for 'Ras_nfpiRaf1'
                                      0.0;		% yinit(13) is the initial condition for 'Ras_GTP'
                                      0.0;		% yinit(14) is the initial condition for 'nfpiRaf1'
                                      Raf1_init_molecules_um_3_mutant;		% yinit(15) is the initial condition for 'Raf1'
                                      0.0;		% yinit(16) is the initial condition for 'Ras_BRaf_pRaf1_tetramer'
                                      0.0;		% yinit(17) is the initial condition for 'Ras_pRaf1_Raf1_tetramer'
                                      0.0;		% yinit(18) is the initial condition for 'nfpRaf1'
                                      0.0;		% yinit(19) is the initial condition for 'Ras_nfpBRaf'
                                      0.0;		% yinit(20) is the initial condition for 'iBRaf'
                                      0.0;		% yinit(21) is the initial condition for 'pMEK'
                                      mE_init_molecules_um_2_mutant;		% yinit(22) is the initial condition for 'mE'
                                      0.0;		% yinit(23) is the initial condition for 'mEL'
                                      0.0;		% yinit(24) is the initial condition for 'pRaf1'
                                      0.0;		% yinit(25) is the initial condition for 'EG2'
                                      0.0;		% yinit(26) is the initial condition for 'Ras_iRaf1'
                                      0.0;		% yinit(27) is the initial condition for 'Ras_nfpiBRaf'
                                      MEK_init_molecules_um_3_mutant;		% yinit(28) is the initial condition for 'MEK'
                                      0.0;		% yinit(29) is the initial condition for 'Ras_Braf_iRaf1_tetramer'
                                      0.0;		% yinit(30) is the initial condition for 'Ras_Raf1'
                                      0.0;		% yinit(31) is the initial condition for 'mELmEL'
                                      0.0;		% yinit(32) is the initial condition for 'nfpiBRaf'
                                      BRaf_init_molecules_um_3_mutant;		% yinit(33) is the initial condition for 'BRaf'
                                      SOS_init_molecules_um_3_mutant;		% yinit(34) is the initial condition for 'SOS'
                                      0.0;		% yinit(35) is the initial condition for 'nfpBRaf'
                                      0.0;		% yinit(36) is the initial condition for 'GRB2_SOS'
                                      GRB2_init_molecules_um_3_mutant;		% yinit(37) is the initial condition for 'GRB2'
                                      0.0;		% yinit(38) is the initial condition for 'Ras_iRaf1_tetramer'
                                      0.0;		% yinit(39) is the initial condition for 'iBRaf_dimer'
                                      Ras_GDP_init_molecules_um_2_mutant;		% yinit(40) is the initial condition for 'Ras_GDP'
                                      0.0;		% yinit(41) is the initial condition for 'E'
                                      0.0;		% yinit(42) is the initial condition for 'Ras_BRaf_Raf1_tetramer'
                                      0.0;		% yinit(43) is the initial condition for 'BRaf_dimer'
                                      0.0;		% yinit(44) is the initial condition for 'EG2SOS'
                                      0.0;		% yinit(45) is the initial condition for 'pERK'
                                      0.0;		% yinit(46) is the initial condition for 'Ras_iRaf1_iBRaf1_tetramer'
];