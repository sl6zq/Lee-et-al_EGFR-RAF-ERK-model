## Define reaction network for extended EGFR-Ras-Raf model:
# Note that cell membrane-bound species are denoted with a subscript m (ₘ)
# to distinguish them from cytoplasmic/volumetric species (no subscript).
using Catalyst, Logging, StructuralIdentifiability
using StatsBase
using PlotlyJS
using DataFrames, DataFramesMeta

#using StructuralIdentifiability, ModelingToolkit
#using ModelingToolkit: t_nounits as t, D_nounits as D



receptor_rn = @reaction_network begin
    # EGF binding (MEMBRANE):
    kEf, EGF + EGFR → mELₘ + EGF    # forward -- EGF assumed to be in excess (and is thus not consumed)
    kEr, EGFR ← mELₘ                # reverse

    # EGFR dimerization (MEMBRANE):
    (kdEf, kdEr), 2mELₘ ↔ mELmELₘ

    # EGFR phosphorylation and dephosphorylation (MEMBRANE):
    (kcatE, kdp), mELmELₘ ↔ Eₘ                  # phosphorylation (forward) and dephosphorylation (reverse)
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
    @observables begin
        totalRASGTP ~ RAStₘ + RAF1ₘ + pRAF1ₘ + nfpRAF1ₘ + BRAFₘ + nfpBRAFₘ
        membRAF ~ RAF1ₘ + pRAF1ₘ + RAF1BRAFₘ + pRAF1BRAFₘ
    end
    # Ras guanine nucleotide exchange rxns (MEMBRANE):
    kRgne/(Kmgne + RAS),  EG2SOSₘ + RAS → RAStₘ + EG2SOSₘ       # GDP-to-GTP exchange (SOS-catalyzed; Michaelis-Menten kinetics)
    kRhydro, RAStₘ → RAS                                        # GTP hydrolysis on Ras

    # Ras-RAF1 binding, all (MEMBRANE):
    (kR1f, kR1r), RAStₘ + RAF1 ↔ RAF1ₘ          # Ras-RAF1, forward & reverse
    (kR1r, kR1f), pRAF1ₘ ↔ RAStₘ + pRAF1        # Ras-pRAF1, reverse & forward
    knfpR1r, nfpRAF1ₘ → RAStₘ + nfpRAF1         # Ras-nfpRAF1, unbinding only

    # Ras-BRAF binding, all (MEMBRANE):
    (kBf, kBr), RAStₘ + BRAF ↔ BRAFₘ            # Ras-BRAF, forward & reverse
    knfpBr, nfpBRAFₘ → RAStₘ + nfpBRAF          # Ras-nfpBRAF, unbinding only

    # RAF1 homodimerization (MEMBRANE):
    (kdR1f, kdR1r), RAF1ₘ + pRAF1ₘ ↔ RAF1pRAF1ₘ         # RAF1-pRAF1 , forward & reverse
    (kdR1r, kdR1f), pRAF1dimₘ ↔ 2pRAF1ₘ                 # pRAF1-pRAF1, reverse & forward

    # BRAF homodimerization (MEMBRANE):
    (kdR1f, kdR1r), 2BRAFₘ ↔ BRAFdimₘ           # BRAF-BRAF, forward & reverse

    # RAF1-BRAF heterodimerization (MEMBRANE):
    (kdR1f, kdR1r), RAF1ₘ + BRAFₘ ↔ RAF1BRAFₘ           # RAF1-BRAF, forward & reverse
    (kdR1r, kdR1f), pRAF1BRAFₘ ↔ pRAF1ₘ + BRAFₘ         # pRAF1-BRAF dimerization, reverse & forward

    # RAF1 phosphorylation (MEMBRANE):
    kpR1, RAF1BRAFₘ → pRAF1BRAFₘ                # BRAF-catalyzed
    kpR1, RAF1pRAF1ₘ → pRAF1dimₘ                # pRAF1-catalyzed

    # pRAF1 dephosphorylation (MEMBRANE and CYTOPLASM):
    kdpR1, pRAF1ₘ → RAF1ₘ                       # pRAF1 dephosphorylation on MEMBRANE
    kdpR1, pRAF1 → RAF1                         # pRAF1 dephosphorylation in CYTOPLASM

    # MEK and ERK dephosphorylation
    kdpMEK, pMEK → MEK
    kdpERK, pERK → ERK

    # SOS dephosphorylation
    kdpSOS, nfpSOS → SOS

    # MEK phosphorylation, catalyzed by active RAF species
    2*kpMEK, BRAFdimₘ + MEK → pMEK + BRAFdimₘ           # BRAFdimₘ-catalzyed
    kpMEK, RAF1BRAFₘ + MEK → pMEK + RAF1BRAFₘ           # RAF1BRAFₘ-catalzyed
    2*kpMEK, pRAF1BRAFₘ + MEK → pMEK + pRAF1BRAFₘ       # pRAF1BRAFₘ-catalzyed
    kpMEK, RAF1pRAF1ₘ + MEK → pMEK + RAF1pRAF1ₘ         # RAF1pRAF1ₘ-catalzyed
    2*kpMEK, pRAF1dimₘ + MEK → pMEK + pRAF1dimₘ         # pRAF1dimₘ-catalzyed
    kpMEK, pRAF1ₘ + MEK → pMEK + pRAF1ₘ                 # pRAF1ₘ-catalzyed
    kpMEK, pRAF1 + MEK → pMEK + pRAF1                   # pRAF1-catalzyed

    # ERK phosphorylation, catalyed by pMEK
    kpERK, pMEK + ERK → pERK + pMEK

    # ERK-mediated negative feedback phosphorylation of membrane-bound species (MEMBRANE):
    knfpBR1, pERK + RAF1ₘ → nfpRAF1ₘ + pERK             # RAF1
    knfpBR1, pERK + BRAFₘ → nfpBRAFₘ + pERK             # BRAF
    knfpSOS, pERK + EG2SOSₘ → EG2ₘ + nfpSOS + pERK      # SOS

end

#obs_def = @reaction_network begin

#end

#"RAF1pRAF1ₘ","pRAF1dimₘ","BRAFdimₘ","RAF1BRAFₘ","pRAF1BRAFₘ"] .* "(t)"

#=
To get species names in the correct order of the original model when
it was run on Julia v1.6.2 (on Rivanna), we have to build the reaction
network backwards, i.e., starting with the last network and adding on
the first network last (as the networks are defined in this script.)
=#
rn = extend(mapk_rn, adaptor_rn)
rn = extend(rn, receptor_rn)
#rn = extend(rn, obs_def)
#rn = extend(rn, var_int)


## Variables for getting model outputs of interest:
pEGFRnames = ["Eₘ","EG2ₘ","EG2SOSₘ"] .* "(t)"
rasgtpnames = ["RAStₘ","RAF1ₘ","pRAF1ₘ","nfpRAF1ₘ","BRAFₘ","nfpBRAFₘ"] .* "(t)"
rasgtpdimnames = ["RAF1pRAF1ₘ","pRAF1dimₘ","BRAFdimₘ","RAF1BRAFₘ","pRAF1BRAFₘ"] .* "(t)"
raf1m_names = ["RAF1ₘ","pRAF1ₘ","RAF1BRAFₘ","pRAF1BRAFₘ"] .* "(t)";  # membrane species containing just one RAF1 molecule
raf1mdim_names = ["RAF1pRAF1ₘ","pRAF1dimₘ"] .* "(t)";    # RAF1 dimers at the membrane (2 RAF1 molecules)

function get_fitting_vars(df; pEGFRnames=pEGFRnames, rasgtpnames=rasgtpnames, rasgtpdimnames=rasgtpdimnames, raf1m_names=raf1m_names, raf1mdim_names=raf1mdim_names)
    rasgtp = sum.(eachrow(df[:, rasgtpnames])) .+ 2.0 .* sum.(eachrow(df[:, rasgtpdimnames]))
    memraf1 = sum.(eachrow(df[:, raf1m_names])) .+ 2.0 .* sum.(eachrow(df[:, raf1mdim_names]))
    pmek = df[:, "pMEK(t)"]
    perk = df[:, "pERK(t)"]

    return rasgtp, memraf1, pmek, perk
end
# convert to ODE system
osys = convert(ODESystem,rn)

# Create equation for observables
@variables y


# which species concentrations can we measure? RAS-GTP, membrane RAF1, pMEK, pERK
obs       = [28, 30]
#count_all = zeros(length(obs))
count_all = DataFrame()
for i in eachindex(obs)
    eq = y ~ states(osys)[obs[i]] #+ states(osys)[30]  # 30:pERK, 28: pMEK
    # are any of the parameters able to be measured?
    # eq = y ~ parameters(osys)[1]
    # Local & global identifiability analysis
    #iden_result = assess_identifiability(osys; measured_quantities=[eq])
    iden_result = assess_local_identifiability(osys; measured_quantities=[eq]) # default p = 0.99
    # extract only the rate constants
    keys_as_vector = collect(keys(iden_result))[33: 64]
    vals_as_vector = collect(values(iden_result))[33: 64]
    # Count occurrences of true/false
    counts_i = countmap(vals_as_vector)
    df = DataFrame(key = collect(keys(counts_i)), count = collect(values(counts_i)), observable = string(states(osys)[obs[i]]) )

    #count_all[i] = counts_i
    append!(count_all, df)
     
end
#count_all[!, :observable] = convert.(String, count_all[:, :observable])

# Create the bar plot
#plot(["Non-identifiable", "Identifiable"], collect(values(count_all)), Layout(title="Observables = pERK",  yaxis_title_text="Number of parameters", uniformtext_minsize=8), kind = "bar")
plot(count_all, x=["Non-identifiable", "Identifiable"], y=:count, color=:observable, kind="bar", Layout(title="Local identifiability", yaxis_title_text="Number of parameters", uniformtext_minsize=7))