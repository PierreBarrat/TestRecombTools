"""
	branch_similarity(t::Tree, MCC1, MCC2)

Return fraction of branches in `t` that are predicted to be in an mcc for both `MCC1` and
  `MCC2`.
"""
function branch_similarity(t::Tree, MCC1, MCC2)
	s = 0
	for n in nodes(t)
		if RecombTools.is_branch_in_mccs(n, MCC1) == RecombTools.is_branch_in_mccs(n, MCC2)
			s += 1
		end
	end
	return s / length(nodes(t))
end


function measure_reproducibility(
	sim, arg::ARG, nit, oa::OptArgs;
	cutoff=0, preresolve=false,
)
	ts = ARGTools.trees_from_ARG(arg)
	trees = Dict(i => t for (i,t) in enumerate(ts))
	if cutoff > 0
		for t in values(trees)
			remove_branches!(t, Exponential(cutoff))
		end
	end
	#
	t1 = trees[1]
	Md = length(t1.lleaves) / nit
	oac = @set oa.Md = Md
	#
	rMCC = ARGTools.MCCs_from_arg(arg)
	MCCs1 = computeMCCs(trees, oac; preresolve)[1,2]
	MCCs2 = computeMCCs(trees, oac; preresolve)[1,2]
	MCCs_naive = computeMCCs(trees, oac; preresolve, naive=true)[1,2]
	return sim(t1, MCCs1, MCCs2), sim(t1, MCCs1, MCCs_naive), sim(t1, MCCs1, rMCC)
end

function eval_reproducibility(nit, N, n, ρ, oa; Nrep = 10, cutoff = 0.)
	dat = zeros(Nrep, 3)
	for rep in 1:Nrep
		arg = ARGTools.SimulateARG.simulate(N, get_r(ρ, n, N, :yule), n)
		dat[rep,:] .= measure_reproducibility(
			branch_similarity, arg, nit, oa;
			cutoff = cutoff*N
		)
	end
	vec(mean(dat, dims=1))
end
