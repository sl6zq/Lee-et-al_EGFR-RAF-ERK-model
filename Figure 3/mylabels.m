names = {'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdE,f'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'};
PLSR_cat = categorical(names);

PLSR_cat = reordercats(PLSR_cat,{'kEf'; 'kcatE'; 'kpMEK'; 'knfpSOS'; 'kBr'; 'kBf'; 'kdpERK'; 'kR1r'; 'kR1f'; 'knfpBR1'; 'kpR1'; 'kSOSr'; 'kSOSf'; 'kRgneslow'; 'kSon'; 'kG2SOSr'; 'kG2SOSf'; 'kiR1r'; 'kiR1f'; 'kSoff'; 'kpERK'; 'kRhydro';
'Kmgneslow'; 'kdpMEK'; 'knfpiR1r'; 'kdR1r'; 'kiBr'; 'kdR1f'; 'kiBf'; 'kfpBr'; 'kdpSOS';
'kdp'; 'kdEr'; 'kdE,f'; 'kdpR1'; 'kG2r'; 'kG2f'; 'knfpiBr'; 'kfpR1r'; 'kEr'});

names_vec = ["kEf" "kcatE" "kpMEK" "knfpSOS" "kBr" "kBf" "kdpERK" "kR1r" "kR1f" "knfpBR1" "kpR1" "kSOSr" "kSOSf" "kRgneslow" "kSon" "kG2SOSr" "kG2SOSf" "kiR1r" "kiR1f" "kSoff" "kpERK" "kRhydro" "Kmgneslow" "kdpMEK" "knfpiR1r" "kdR1r" "kiBr" "kdR1f" "kiBf" "kfpBr" "kdpSOS" "kdp" "kdEr" "kdE,f" "kdpR1" "kG2r" "kG2f" "knfpiBr" "kfpR1r" "kEr"]';
