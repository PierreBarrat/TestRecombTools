"""
    consistent_mcc_triplets(M12, M13, M23)
"""
function consistent_mcc_triplets(M12, M13, M23, trees)
	s = 0
	Z = 0
	for n in nodes(trees[1])
		if RecombTools.is_branch_in_mccs(n, M12) && RecombTools.is_branch_in_mccs(n, M13)
			Z += 1
			# Must find a corresponding common branch in trees[2] (or trees[3])
			# We find it by taking the intesection of the clade of n with the MCC n
			# belongs to
			m12 = RecombTools.find_mcc_with_branch(n, M12)[2]
			m13 = RecombTools.find_mcc_with_branch(n, M13)[2]
			clade = [x.label for x in POTleaves(n)]
			n2 = lca(trees[2], intersect(clade, m12))
			n3 = lca(trees[3], intersect(clade, m13))
			if !RecombTools.is_branch_in_mccs(n2, M23) || !RecombTools.is_branch_in_mccs(n3, M23)
				s += 1
			end
		end
	end
	return Z == 0 ? 0.0 : s / Z
end



#####################
##################### REMOVING BRANCHES
#####################
"""
    remove_branches!(n::TreeNode, p::Distribution)
"""
function remove_branches!(n::TreeNode, p::Distribution)
    child_list = copy(n.child)
	r = rand(p)
    if !ismissing(n.tau) && n.tau < r
        if !n.isleaf
            nr = delete_node!(n)
        else
            n.tau = 0.
        end
    end

    for c in child_list
        remove_branches!(c, p)
    end

    return nothing
end

"""
    remove_branches!(t::Tree, p)

Stochastically remove branches from `t` using probability distribution `p`
"""
function remove_branches!(t::Tree, p)
    remove_branches!(t.root, p)
    node2tree!(t, t.root)
    nothing
end
function remove_branches!(t::Tree, c::Real, N::Int)
	if c > 0
		remove_branches!(t, Exponential(c*N))
	end
	return nothing
end

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

"""
	varinfo_similarity(MCCs...)
"""
varinfo_similarity(t::Tree, MCCs...) = varinfo_similarity(MCCs...)
function varinfo_similarity(MCCs...)
	leaves = sort(vcat(first(MCCs)...))
	assignments = [assignment_vector(leaves, mccs) for mccs in MCCs]
	out = 0
	Z = 0
	for i in 1:length(MCCs), j in (i+1):length(MCCs)
		out += Clustering.varinfo(assignments[i], assignments[j])
		Z += 1
	end
	return out / Z
end
"""
"""
function rand_index_similarity(MCCs...)
	leaves = sort(vcat(first(MCCs)...))
	assignments = [assignment_vector(leaves, mccs) for mccs in MCCs]
	out = 0
	Z = 0
	for i in 1:length(MCCs), j in (i+1):length(MCCs)
		out += (1 - Clustering.randindex(assignments[i], assignments[j])[2])
		Z += 1
	end
	return out / Z
end
function assignment_vector(leaves, MCCs)
	a = zeros(Int, length(leaves))
	for (k,m) in enumerate(MCCs)
		for n in m
			a[findfirst(==(n), leaves)] = k
		end
	end

	return a
end

function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / N
    elseif simtype == :yule
        return ρ / N
    else
        @error "Unrecognized `simtype`."
    end
end
