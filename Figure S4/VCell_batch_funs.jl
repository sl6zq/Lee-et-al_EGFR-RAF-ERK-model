## Packages:
using DelimitedFiles

## Functions:
"""
Load the results of batch simulation generated from VCell.\n
-------\n
- `params` = Matrix or DataFrame of parameter values with parameters as columns and parameter sets as rows.
- `dir` = The directory/path containing the batch simulation results from VCell.
"""
function loadbatch(params, dir_batchres)

    n = size(params, 1)    # total number of batch simulations
    batchres = []   # initialize empty array for storing simulation results
    inds = []   # array/vector for keeping track of which simulations had completed results
    vars = []   # initialize empty array for names of model variables (species)
    # Load the batch results:
    for i in 1:n
        j = i - 1
        # Define name of file containing results from simulation i:
        if j < 10
            sim_name = "00"*string(j)
        elseif j < 100
            sim_name = "0"*string(j)
        else
            sim_name = string(j)
        end
        sim_path = dir_batchres * "/" * sim_name * ".txt"
        if isfile(sim_path)
            res, var = readdlm(sim_path, header=true)
            if isempty(vars)
                vars = vec(var)   # store variable names (as a vector)
            end
            if size(res, 1) < 2
                break
            end
            push!(batchres, res)    # store simulation results
            push!(inds, i)
        end
    end
    return batchres, vars, inds
end


function loadbatch2(params, dir_batchres, n_reps)   # for loading stochasticity results with repeats from raw text files provided by SL
    p = size(params, 1)    # total number of batch simulations
    batchres = []   # initialize empty array for storing simulation results
    # Load the batch results:
    for i in 1:p
        j = i - 1   # index for keeping track of numbering in file names
        batchres_i = []
        if j < 100
            sim_name = "0"*string(j)*" "
        else
            sim_name = string(j)
        end
        for k in 1:n_reps
            sim_path = dir_batchres * "/" * sim_name * string(k) * ".txt"   # name of file containing results from simulation i and replicate k
            if isfile(sim_path)
                res_k = CSV.read(sim_path, DataFrame)                       # load VCell simulation results
                rename!(res_k, [replace(names(res_k)[i], r"_Count" => "") for i in eachindex(names(res_k))])    # cleaning up species names
                push!(batchres_i, res_k)                                    # store VCell simulation results
            end
        end
        push!(batchres, batchres_i)
    end
    return batchres
end

## Dictionaries for model outputs of interest:
# Dictionary for number of pRAF1 molecules per species:
pRAF1_dict = Dict{String, Integer}(
    "pRaf1" => 1,
    "Ras_pRaf1_Raf1_tetramer" => 1,
    "Ras_pRaf1" => 1,
    "Ras_pRaf1_tetramer" => 2,
    "Ras_BRaf_pRaf1_tetramer" => 1,
    "Ras_pRaf1_iBRaf_tetramer" => 1
    )

# Dictionary for number of RAF1 molecules on membrane-bound species:
memRAF1_dict = Dict{String, Integer}(
    "Ras_Raf1" => 1,
    "Ras_iRaf1" => 1,
    "Ras_pRaf1_Raf1_tetramer" => 2,
    "Ras_pRaf1" => 1,
    "Ras_BRaf_Raf1_tetramer" => 1,
    "Ras_iRaf1_tetramer" => 2,
    "Ras_pRaf1_tetramer" => 2,
    "Ras_BRaf_pRaf1_tetramer" => 1,
    "Ras_pRaf1_iBRaf_tetramer" => 1,
    "Ras_Braf_iRaf1_tetramer" => 1,
    "Ras_Raf1_iBRaf_tetramer" => 1,
    "Ras_nfpRaf1" => 1,
    "Ras_nfpiRaf1" => 1,
    "Ras_iRaf1_iBRaf1_tetramer" => 1
)

# Dictionary for total number of Ras-GTP molecules:
RASt_dict = Dict{String, Integer}(
    "Ras_GTP" => 1,
    "Ras_BRaf" => 1,
    "Ras_iBRaf" => 1,
    "Ras_nfpBRaf" => 1,
    "Ras_nfpiBRaf" => 1,
    "Ras_iRaf1_iBRaf1_tetramer" => 2, ################################## spelling on BRaf!!!
    "Ras_Raf1" => 1,
    "Ras_iRaf1" => 1,
    "Ras_pRaf1_Raf1_tetramer" => 2,
    "Ras_pRaf1" => 1,
    "Ras_BRaf_Raf1_tetramer" => 2,
    "Ras_iRaf1_tetramer" => 2,
    "Ras_pRaf1_tetramer" => 2,
    "Ras_BRaf_pRaf1_tetramer" => 2,
    "Ras_pRaf1_iBRaf_tetramer" => 2,
    "Ras_Braf_iRaf1_tetramer" => 2,
    "Ras_Raf1_iBRaf_tetramer" => 2,
    "Ras_nfpRaf1" => 1,
    "Ras_nfpiRaf1" => 1,
    "Ras_iRaf1_iBRaf1_tetramer" => 2
)

# Dictionary for total phosphorylated EGFR:
pEGFR_dict = Dict{String, Integer}(
    "E" => 2,
    "EG2" => 2,
    "EG2SOS" => 2
)

# Function for extracting fitting variables of interest:
function get_fitting_vars_vcell(df;
                                pRAF1_dict = pRAF1_dict,
                                memRAF1_dict = memRAF1_dict,
                                RASt_dict = RASt_dict,
                                pEGFR_dict = pEGFR_dict)
    # Membrane-bound RAF1:
    memRAF1_names = string.(keys(memRAF1_dict))   # names of species (columns) containing membrane RAF1
    memRAF1_molec = [get(memRAF1_dict, i, 0) for i in memRAF1_names] |> permutedims  # number of membrane RAF1 molecules per species
    memraf1 = sum.(eachrow(memRAF1_molec.*df[:, memRAF1_names]));    # max values of total membrane RAF1 per simulation

    # Total Ras-GTP:
    RASt_names = string.(keys(RASt_dict))  # names of species (columns) containing Ras-GTP
    RASt_molec = [get(RASt_dict, i, 0) for i in RASt_names] |> permutedims  # number of Ras-GTP molecules per species
    rasgtp = sum.(eachrow(RASt_molec.*df[:, RASt_names]));    # max values of total Ras-GTP per simulation

    # pMEK and pERK:
    pmek = df[:, "pMEK"]
    perk = df[:, "pERK"]

    return rasgtp, memraf1, pmek, perk
    # return (rasgtp=rasgtp, memraf1=memraf1, pmek=pmek, perk=perk)
end