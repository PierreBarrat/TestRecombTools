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
    if !ismissing(n.tau) && n.tau < rand(p)
        if !n.isleaf
            nr = delete_node!(n)
            for c in nr.child
                remove_branches!(c, p)
            end
        else
            n.tau = 0.
        end
    else
        for c in n.child
            remove_branches!(c, p)
        end
    end
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

