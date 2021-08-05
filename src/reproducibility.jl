


function measure_reproducibility(
	sim, arg::ARG, nit, oa::OptArgs;
	cutoff=0, preresolve=false, Nrep = 10,
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
	sim = varinfo_similarity
	dat = zeros(Nrep, 3)
	for rep in 1:Nrep
		arg = ARGTools.SimulateARG.simulate(N, get_r(ρ, n, N, :yule), n)
		dat[rep,:] .= measure_reproducibility(
			sim, arg, nit, oa;
			cutoff = cutoff*N
		)
	end
	vec(mean(dat, dims=1))
end

"""
	remove_low_bootstrap!(tree::Tree, bmin)
"""
function remove_low_bootstrap!(tree::Tree, bmin)
	for n in internals(tree)
		if !n.isroot
			bs = tryparse(Int, split(n.label, "__")[1])
			if !isnothing(bs) && bs < bmin
				delete_node!(n)
			end
		end
	end

	node2tree!(tree, tree.root)
end


function remove_low_bootstrap(t::Tree, bmin)
	tree = copy(t)
	remove_low_bootstrap!(tree, bmin)
	return tree
end
function remove_low_bootstrap(treefile::AbstractString, bmin)
	t = read_tree(treefile)
	remove_low_bootstrap!(t, bmin)
	return t
end
