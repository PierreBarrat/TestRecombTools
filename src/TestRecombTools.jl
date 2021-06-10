module TestRecombTools

using ARGTools
using Clustering
using CSV
using DataFrames
using Distributions
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

include("reproducibility.jl")
export eval_reproducibility

include("eval_gamma.jl")
export eval_gamma

end # module
