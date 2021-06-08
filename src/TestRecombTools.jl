module TestRecombTools

using ARGTools
using DataFrames
using Distributions
using RecombTools
using StatsBase
using TreeTools

export TRT
const TRT = TestRecombTools

include("main.jl")
export eval_mcc_inf
export eval_runopt
export remove_branches!

end # module
