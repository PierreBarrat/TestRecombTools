module TestRecombTools

using ARGTools
using Clustering
using CSV
using DataFrames
using Distributions
using JSON3
using RecombTools
using Setfield
using StatsBase
using TreeTools

export TRT
const TRT = TestRecombTools

include("main.jl")
export eval_mcc_inf
export eval_runopt
export remove_branches!

include("tools.jl")
export consistent_mcc_triplets

include("reproducibility.jl")
export eval_reproducibility

include("eval_gamma.jl")
export eval_gamma

include("clustering_distance.jl")
include("splits.jl")
include("eval_branches.jl")
include("reassortments.jl")
include("interaction_w_auspice.jl")

end # module
