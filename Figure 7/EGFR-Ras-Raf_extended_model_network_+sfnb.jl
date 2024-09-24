## Define reaction network for extended EGFR-Ras-Raf model:
# Note that cell membrane-bound species are denoted with a subscript m (ₘ)
# to distinguish them from cytoplasmic/volumetric species (no subscript).
using Catalyst
receptor_rn = @reaction_network begin
    # EGF binding (MEMBRANE):
    kEf, EGF + EGFR → mELₘ + EGF        # forward -- EGF assumed to be in excess (and is thus not consumed)
    kEr, EGFR ← mELₘ                    # reverse

    # EGFR dimerization (MEMBRANE):
    (kdEf, kdEr), 2mELₘ ↔ mELmELₘ

    # EGFR phosphorylation and dephosphorylation (MEMBRANE):
    (kcatE, kdp), mELmELₘ ↔ Eₘ          # phosphorylation (forward) and dephosphorylation (reverse)
end

adaptor_rn = @reaction_network begin
    # EGFR-GRB2 binding (MEMBRANE):
    (kG2f, kG2r), Eₘ + GRB2 ↔ EG2ₘ

    # EGFR-GRB2:SOS binding (MEMBRANE):
    (kG2SOSf, kG2SOSr), Eₘ + G2SOS ↔ EG2SOSₘ

    # EGFR:GRB2-SOS binding (MEMBRANE):
    (kSOSf, kSOSr), EG2ₘ + SOS ↔ EG2SOSₘ

    # GRB2-SOS binding
    (kSOSf, kSOSr), GRB2 + SOS ↔ G2SOS
end

mapk_rn = @reaction_network begin
    # Ras guanine nucleotide exchange rxns (MEMBRANE):
    kRgne*EG2SOSₘ/(Kmgne + RAS),  RAS → RAStₘ       # GDP-to-GTP exchange (SOS-catalyzed; Michaelis-Menten kinetics)
    kRhydro, RAStₘ → RAS                            # GTP hydrolysis on Ras

    # Sorafenib binding/unbinding with RAF1/BRAF
    kSon, (RAF1 + sfnb, BRAF + sfnb) → (iRAF1 + sfnb, iBRAF + sfnb)     # sorafenib binding, forward
    kSoff, (iRAF1, iBRAF) → (RAF1, BRAF)                                # sorafenib binding, reverse (i.e., unbinding)

    # Ras-RAF1 binding, all (MEMBRANE):
    (kR1f, kR1r), RAStₘ + RAF1 ↔ RAF1ₘ          # Ras-RAF1, forward & reverse
    (kiR1f, kiR1r), RAStₘ + iRAF1 ↔ iRAF1ₘ      # Ras-iRAF1, forward & reverse
    (kR1r, kR1f), pRAF1ₘ ↔ RAStₘ + pRAF1        # Ras-pRAF1, reverse & forward
    knfpR1r, nfpRAF1ₘ → RAStₘ + nfpRAF1         # Ras-nfpRAF1, unbinding only
    knfpiR1r, nfpiRAF1ₘ → RAStₘ + nfpiRAF1      # Ras-nfpiRAF1, unbinding only

    # Ras-BRAF binding, all (MEMBRANE):
    (kBf, kBr), RAStₘ + BRAF ↔ BRAFₘ            # Ras-BRAF, forward & reverse
    (kiBf, kiBr), RAStₘ + iBRAF ↔ iBRAFₘ        # Ras-iBRAF, forward & reverse
    knfpBr, nfpBRAFₘ → RAStₘ + nfpBRAF          # Ras-nfpBRAF, unbinding only
    knfpiBr, nfpiBRAFₘ → RAStₘ + nfpiBRAF       # Ras-nfpiBRAF, unbinding only

    # RAF1 homodimerization (MEMBRANE):
    (kdR1f, kdR1r), RAF1ₘ + pRAF1ₘ ↔ RAF1pRAF1ₘ     # RAF1-pRAF1 , forward & reverse
    (kdR1r, kdR1f), pRAF1dimₘ ↔ 2pRAF1ₘ             # pRAF1-pRAF1, reverse & forward
    (kdR1f, kdR1r), 2iRAF1ₘ ↔ iRAF1dimₘ             # iRAF1-iRAF1, forward & reverse

    # BRAF homodimerization (MEMBRANE):
    (kdR1f, kdR1r), 2BRAFₘ ↔ BRAFdimₘ               # BRAF-BRAF, forward & reverse
    (kdR1f, kdR1r), BRAFₘ + iBRAFₘ ↔ BRAFiBRAFₘ     # BRAF-iBRAF, forward & reverse
    (kdR1f, kdR1r), 2iBRAFₘ ↔ iBRAFdimₘ             # iBRAF-iBRAF, forward & reverse

    # RAF1-BRAF heterodimerization (MEMBRANE):
    (kdR1f, kdR1r), RAF1ₘ + BRAFₘ ↔ RAF1BRAFₘ           # RAF1-BRAF, forward & reverse
    (kdR1f, kdR1r), RAF1ₘ + iBRAFₘ ↔ RAF1iBRAFₘ         # RAF1-iBRAF, forward & reverse
    (kdR1f, kdR1r), iRAF1ₘ + BRAFₘ ↔ iRAF1BRAFₘ         # iRAF-BRAF, forward & reverse
    (kdR1f, kdR1r), iRAF1ₘ + iBRAFₘ ↔ iRAF1iBRAFₘ       # iRAF1-iBRAF, forward & reverse
    (kdR1r, kdR1f), pRAF1BRAFₘ ↔ pRAF1ₘ + BRAFₘ         # pRAF1-BRAF dimerization, reverse & forward
    (kdR1r, kdR1f), pRAF1iBRAFₘ ↔ pRAF1ₘ + iBRAFₘ       # pRAF1-iBRAF dimerization, reverse & forward

    # RAF1 phosphorylation (MEMBRANE):
    kpR1, RAF1BRAFₘ → pRAF1BRAFₘ      # BRAF-catalyzed
    kpR1, RAF1pRAF1ₘ → pRAF1dimₘ      # pRAF1-catalyzed
    kpR1, RAF1iBRAFₘ → pRAF1iBRAFₘ    # iBRAF-catalyzed

    # pRAF1 dephosphorylation (MEMBRANE and CYTOPLASM):
    kdpR1, pRAF1ₘ → RAF1ₘ       # pRAF1 dephosphorylation on MEMBRANE
    kdpR1, pRAF1 → RAF1         # pRAF1 dephosphorylation in CYTOPLASM

    # MEK and ERK dephosphorylation
    kdpMEK, pMEK → MEK
    kdpERK, pERK → ERK

    # SOS dephosphorylation
    kdpSOS, nfpSOS → SOS

    # MEK phosphorylation, catalyzed by active RAF species (mass-action)
    2*kpMEK, BRAFdimₘ + MEK → pMEK + BRAFdimₘ           # BRAFdimₘ-catalzyed
    kpMEK, BRAFiBRAFₘ + MEK → pMEK + BRAFiBRAFₘ         # BRAFiBRAFₘ-catalzyed
    kpMEK, RAF1BRAFₘ + MEK → pMEK + RAF1BRAFₘ           # RAF1BRAFₘ-catalzyed
    kpMEK, iRAF1BRAFₘ + MEK → pMEK + iRAF1BRAFₘ         # iRAF1BRAFₘ-catalzyed
    2*kpMEK, pRAF1BRAFₘ + MEK → pMEK + pRAF1BRAFₘ       # pRAF1BRAFₘ-catalzyed
    kpMEK, pRAF1iBRAFₘ + MEK → pMEK + pRAF1iBRAFₘ       # pRAF1iBRAFₘ-catalzyed
    kpMEK, RAF1pRAF1ₘ + MEK → pMEK + RAF1pRAF1ₘ         # RAF1pRAF1ₘ-catalzyed
    2*kpMEK, pRAF1dimₘ + MEK → pMEK + pRAF1dimₘ         # pRAF1dimₘ-catalzyed
    kpMEK, pRAF1ₘ + MEK → pMEK + pRAF1ₘ                 # pRAF1ₘ-catalzyed
    kpMEK, pRAF1 + MEK → pMEK + pRAF1                   # pRAF1-catalzyed

    # ERK phosphorylation, catalyed by pMEK
    kpERK, pMEK + ERK → pERK + pMEK

    # ERK-mediated negative feedback phosphorylation of membrane-bound species (mass-action) (MEMBRANE):
    knfpBR1, pERK + RAF1ₘ → nfpRAF1ₘ + pERK             # RAF1
    knfpBR1, pERK + BRAFₘ → nfpBRAFₘ + pERK             # BRAF
    knfpBR1, pERK + iRAF1ₘ → nfpiRAF1ₘ + pERK           # iRAF1
    knfpBR1, pERK + iBRAFₘ → nfpiBRAFₘ + pERK           # iBRAF
    knfpSOS, pERK + EG2SOSₘ → EG2ₘ + nfpSOS + pERK      # SOS
end

#=
To get species names in the correct order of the original model when
it was run on Julia v1.6.2 (on Rivanna), we have to build the reaction
network backwards, i.e., starting with the last network and adding on
the first network last (as the networks are defined in this script.)
=#
rn = extend(mapk_rn, adaptor_rn)
rn = extend(rn, receptor_rn)

## Variables for getting model outputs of interest:
pEGFRnames = ["Eₘ","EG2ₘ","EG2SOSₘ"] .* "(t)"
rasgtpnames = ["RAStₘ","RAF1ₘ","pRAF1ₘ","nfpRAF1ₘ","nfpiRAF1ₘ","BRAFₘ","nfpBRAFₘ","nfpiBRAFₘ","iRAF1ₘ","iBRAFₘ"] .* "(t)"
rasgtpdimnames = ["RAF1pRAF1ₘ","pRAF1dimₘ","iRAF1dimₘ",
    "BRAFdimₘ","iBRAFdimₘ","BRAFiBRAFₘ",
    "RAF1BRAFₘ","iRAF1BRAFₘ","iRAF1iBRAFₘ","pRAF1BRAFₘ","pRAF1iBRAFₘ","RAF1iBRAFₘ"] .* "(t)"
raf1m_names = ["RAF1ₘ","iRAF1ₘ","pRAF1ₘ","RAF1BRAFₘ","pRAF1BRAFₘ","iRAF1BRAFₘ","iRAF1iBRAFₘ","pRAF1iBRAFₘ","RAF1iBRAFₘ"] .* "(t)";  # membrane species containing just one RAF1 molecule
raf1mdim_names = ["RAF1pRAF1ₘ","pRAF1dimₘ","iRAF1dimₘ"] .* "(t)";    # RAF1 dimers at the membrane (2 RAF1 molecules)

function get_fitting_vars(df; pEGFRnames=pEGFRnames, rasgtpnames=rasgtpnames, rasgtpdimnames=rasgtpdimnames, raf1m_names=raf1m_names, raf1mdim_names=raf1mdim_names)
    rasgtp = sum.(eachrow(df[:,rasgtpnames])) .+ 2.0 .* sum.(eachrow(df[:,rasgtpdimnames]))
    memraf1 = sum.(eachrow(df[:,raf1m_names])) .+ 2.0 .* sum.(eachrow(df[:,raf1mdim_names]))
    pmek = df[:, "pMEK(t)"]
    perk = df[:, "pERK(t)"]

    return rasgtp, memraf1, pmek, perk
end