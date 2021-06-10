function distance_to_naive(
	sim, arg::ARG, oa::OptArgs, γvals;
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
	MCCs_naive = computeMCCs(trees, oa; preresolve, naive=true)[1,2]
	rMCC = ARGTools.MCCs_from_arg(arg)

	dat = Array{Float64, 2}(undef, 0, 3)
	for γ in γvals
		oac = @set oa.γ = γ
		iMCCs = computeMCCs(trees, oac; preresolve)[1,2]
		dat = [dat ; [γ sim(trees[1], iMCCs, MCCs_naive) sim(trees[1], iMCCs, rMCC)]]
	end

	return dat
end


function eval_gamma(γvals, N, n, ρ, oa; Nrep = 10, cutoff = 0.)
	sim = varinfo_similarity
	dat = zeros(Float64, length(γvals), 3)
	for rep in 1:Nrep
		arg = ARGTools.SimulateARG.simulate(N, get_r(ρ, n, N, :yule), n)
		dat .+= distance_to_naive(sim, arg,oa, γvals; cutoff)
	end

	return dat / Nrep
end
