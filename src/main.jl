#####################
##################### Simulating ARGs and computing MCCs. Writing result.
#####################


"""
	simulate(;
		# ARG simulation
		N = 100_000,
		n = 100,
		ρ = 0.1,
		simtype = :yule,
		cutoff = 0.,
		K = 2,
		nmax = n,
		s = (nmax - n) / (2*N),
		# Inference
		γ = 2,
		nit = 100,
		resolve = true,
		# Others
		Nrep = 250,
		write = true,
		verbose = false,
		outfolder = make_outfolder_name(N, n, ρ, cutoff, γ, nit, resolve),
	)
"""
function simulate(;
	# ARG simulation
	N = 100_000,
	n = 100,
	ρ = 0.1,
	simtype = :yule,
	cutoff = 0.,
	K = 2,
	nmax = n,
	s = (nmax - n) / (2*N),
	# Inference
	γ = 2,
	nit = 100,
	resolve = true,
	# Others
	Nrep = 250,
	write = true,
	verbose = false,
	outfolder = make_outfolder_name(N, n, ρ, cutoff, γ, nit, resolve),
)
	for rep in 1:Nrep
		verbose && println("Simulating trees ... ")
		arg, trees = simulate_trees(N, n, ρ, simtype, cutoff, K, nmax, s)
		verbose && println("Inferring MCCs ... ")
		allMCCs = get_mccs(arg, trees, γ, nit, resolve; verbose) # array of length 4
		write && write_mccs(outfolder, rep, allMCCs, "a")
		write && write_trees(outfolder, rep, trees, "a")
	end

	return nothing
end

function write_trees(outfolder, id, trees, mode="a")
	dir_id = 1
	# dir_id = id
	mkpath(outfolder * "/rep$(dir_id)")

	for (i,t) in enumerate(trees)
		i > 2 && break
		of = outfolder * "/rep$(dir_id)/tree$(i).nwk"
		open(of, mode) do w
			id > 1 ? write(w, "\n# rep $(id)\n") : write(w, "# rep $(id)\n")
		end
		TreeTools.write_newick(of, t, "a")
	end
end

function write_mccs(outfolder, id, allMCCs, mode="a")
	dir_id = 1
	# dir_id = id
	mkpath(outfolder * "/rep$(dir_id)")

	# Real
	of = outfolder * "/rep$(dir_id)/realMCCs.dat"
	open(of, mode) do w
		id > 1 ? write(w, "\n# rep $(id)\n") : write(w, "# rep $(id)\n")
	end
	TreeKnit.write_mccs(of, allMCCs[1], "a")

	# Naive
	of = outfolder * "/rep$(dir_id)/naiveMCCs.dat"
	open(of, mode) do w
		id > 1 ? write(w, "\n# rep $(id)\n") : write(w, "# rep $(id)\n")
	end
	TreeKnit.write_mccs(of, allMCCs[2], "a")

	# Inferred
	of = outfolder * "/rep$(dir_id)/inferredMCCs.dat"
	open(of, mode) do w
		id > 1 ? write(w, "\n# rep $(id)\n") : write(w, "# rep $(id)\n")
	end
	TreeKnit.write_mccs(of, allMCCs[3], "a")
end

function get_mccs(arg, trees, γ, nit, resolve; verbose)
	oa = OptArgs(; γ, resolve, nMCMC = nit, verbose)

	rMCCs = ARGTools.MCCs_from_arg(arg, 1, 2)
	iMCCs1 = computeMCCs(trees[1], trees[2], oa)
	naiveMCCs = computeMCCs(trees[1], trees[2], oa; naive=true)

	return rMCCs, naiveMCCs, iMCCs1 #, iMCCs2
end

function simulate_trees(N, n, ρ, simtype, cutoff, K, nmax, s)
	# ARG
	arg = ARGTools.SimulateARG.simulate(N, get_r(ρ, n, N, simtype), n; K, simtype, nmax, s)
	# Trees
	trees = ARGTools.trees_from_ARG(arg)
	# trees = Dict(i => t for (i,t) in enumerate(ts))
	if cutoff > 0
		for t in trees
			remove_branches!(t, Exponential(cutoff * N))
		end
	end

	return arg, trees
end







#####################
##################### Evaluating MCC inference
#####################



"""
    eval_naive_inf(N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
        cutoff = 0.,
        Nrep = 1,
        sfields::Tuple = (:ρ,:cutoff),
        out = ""
    )
"""
function eval_naive_inf(N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    cutoff = 0.,
    Nrep = 1,
    sfields::Tuple = (:ρ,:cutoff),
    out = ""
)
    #
    args = Dict(:N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype, :cutoff=>cutoff)
    #
    dat = DataFrame(df_fields())
    for f in sfields
        insertcols!(dat, f=>Any[])
    end
    #
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        t1, t2 = ARGTools.trees_from_ARG(arg)
        if cutoff != 0
            remove_branches!(t1, Distributions.Exponential(cutoff*N))
            remove_branches!(t2, Distributions.Exponential(cutoff*N))
        end
        trees = Dict(1=>t1, 2=>t2)
        MCCs = computeMCCs!(trees; naive=true)
        return MCCs
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        args[:time] = td[2]
        args[:bytes] = td[3]
        for f in sfields
            d[f] = args[f]
        end
        push!(dat, d)
    end
    #
    out != "" && CSV.write(out, dat)
    return dat
end

"""
    eval_real(N::Int, n::Int, ρ::Float64, simtype::Symbol;
        Nrep = 1,
        out = "",
    )
"""
function eval_real(N::Int, n::Int, ρ::Float64, simtype::Symbol;
    Nrep = 1,
    out = "",
)
    dat = DataFrame(df_fields())
    insertcols!(dat, :ρ=>Float64[])
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        t1, t2 = ARGTools.trees_from_ARG(arg)
        MCCs = Dict((1,2) => ARGTools.MCCs_from_arg(arg))
        return MCCs
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        d[:ρ] = ρ
        push!(dat, d)
    end
    #
    out != "" && CSV.write(out, dat)
    return dat
end


"""
    eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    Md=10, lk_sort=true,
    cutoff = 0.,
    Nrep = 1,
    sfields::Tuple = (:ρ,:cutoff)
    out = "",
    verbose=false
)

Eval the performance of `SplitGraph.runopt` at inferring MCCs.
"""
function eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    Md=10, lk_sort=true,
    cutoff = 0.,
    Nrep = 1,
    sfields::Tuple = (:ρ,:cutoff),
    out = "",
    verbose=false
)
    #
    args = Dict(:γ=>γ, :N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype,
        :Md=>Md, :lk_sort=>lk_sort,
        :cutoff=>cutoff,
    )
    #
    dat = DataFrame(df_fields())
    for f in sfields
        insertcols!(dat, f=>Any[])
    end
    #
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        let cutoff=cutoff, N=N
            t1, t2 = ARGTools.trees_from_ARG(arg)
            trees = Dict(1=>t1, 2=>t2)
            if cutoff != 0
                remove_branches!(t1, Distributions.Exponential(cutoff*N))
                remove_branches!(t2, Distributions.Exponential(cutoff*N))
            end
            oa = OptArgs(
                γ=γ,
                Md=Md,
                likelihood_sort=lk_sort,
                verbose=verbose,
            )
            try
                MCCs = computeMCCs!(
                    trees, oa;
                )
                return MCCs
            catch err
                mkpath("tmp")
                write_newick("tmp/tree1.nwk", trees[1])
                write_newick("tmp/tree2.nwk", trees[2])
                error(err)
            end
        end
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        args[:time] = td[2]
        args[:bytes] = td[3]
        for f in sfields
            d[f] = args[f]
        end
        push!(dat, d)
    end
    #
    out != "" && CSV.write(out, dat)
    return dat
end

"""
    eval_runopt_manytrees(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
        K = 3,
        Md=10, lk_sort=true,
        cutoff = 0.,
        Nrep = 1,
        sfields::Tuple = (:ρ,:cutoff),
        out = "",
        verbose=false
    )

Eval the performance of `SplitGraph.runopt` at inferring MCCs.
"""
function eval_runopt_manytrees(γ::Real, N::Int, n::Int, ρ::Float64, simtype::Symbol, K::Int;
     Md=10, lk_sort=true,
    cutoff = 0.,
    Nrep = 1,
    sfields::Tuple = (:ρ,:cutoff),
    out = "",
    verbose=false
)
    #
    args = Dict(:γ=>γ, :N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype,
        :Md=>Md, :lk_sort=>lk_sort,
        :cutoff=>cutoff,
    )
    #
    dat = DataFrame(df_fields())
    for f in sfields
        insertcols!(dat, f=>Any[])
    end
    #
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        let cutoff=cutoff, N=N
            ts = ARGTools.trees_from_ARG(arg)
            trees = Dict(i=>t for (i,t) in enumerate(ts))
            if cutoff != 0
                for t in values(trees)
                    remove_branches!(t, Distributions.Exponential(cutoff*N))
                end
            end
            oa = OptArgs(
                γ=γ,
                Md=Md,
                likelihood_sort=lk_sort,
                verbose=verbose,
            )
            try
                MCCs = computeMCCs!(trees, oa;)
                return MCCs
            catch err
                mkpath("tmp")
                write_newick("tmp/tree1.nwk", trees[1])
                write_newick("tmp/tree2.nwk", trees[2])
                error(err)
            end
        end
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r; simtype, K)
        d = td[1]
        args[:time] = td[2]
        args[:bytes] = td[3]
        for f in sfields
            d[f] = args[f]
        end
        push!(dat, d)
    end
    #
    out != "" && CSV.write(out, dat)
    return dat
end

function simulate_sequences!(trees::Dict, N, c = 0.5, L = 1000; polytomies=true)
    μ = c == 0 ? 1 / (N*L*0.2) * 1.25 : 1 / (N*L*c) * 1.25
    # Simulate sequences && prune branches with no mutations
    for (i,t) in trees
        TreeAlgs.Evolve.evolve!(t, L, μ;  seqkey = :selfseq)
        TreeTools.compute_mutations!(n -> n.data.dat[:selfseq], t, :realmuts)
        c != 0 && polytomies && TreeTools.delete_branches!(n->isempty(n.data.dat[:realmuts]), t)
    end
    #
    return nothing
end

function _eval_mcc_inf(rMCC, iMCC, t1::Tree)
    # Mean MCC sizes
    rmMCC = mean([length(x) for x in values(rMCC)])
    imMCC = mean([length(x) for x in values(iMCC)])
    # Common branches
    τc_p = 0.; τnc_p = 0.; τc_n = 0.; τnc_n = 0.
    νc_p = 0.; νnc_p = 0.; νc_n = 0.; νnc_n = 0.
    for n in Iterators.filter(x->!x.isroot, values(t1.lnodes))
        if TreeKnit.is_branch_in_mccs(n,iMCC) # Branch predicted to be shared with the other tree
            if TreeKnit.is_branch_in_mccs(n, rMCC) # Correct prediction
                τc_p += n.tau
                νc_p += 1.
            else
                τc_n += n.tau
                νc_n += 1.
            end
        else # Branch predicted to not be shared with the other tree
            if !TreeKnit.is_branch_in_mccs(n, rMCC) # Correct prediction
                τnc_p += n.tau
                νnc_p += 1.
            else
                τnc_n += n.tau
                νnc_n += 1.
            end
        end
    end
    νc_p /= (length(t1.lnodes)-1); νnc_p /= (length(t1.lnodes)-1); νc_n /= (length(t1.lnodes)-1); νnc_n /= (length(t1.lnodes)-1)
    T = sum(skipmissing(x.tau for x in values(t1.lnodes))) # Total branch length
    τc_p /= T; τnc_p /= T; τc_n /= T; τnc_n /= T

    return Dict(:τc_p => τc_p, :τnc_p => τnc_p, :τc_n => τc_n, :τnc_n => τnc_n,
                :νc_p => νc_p, :νnc_p => νnc_p, :νc_n => νc_n, :νnc_n => νnc_n,
                :mMCCsize=>imMCC, :nMCC=>length(iMCC))
end

function _eval_split_inf(true_splits, resolved_splits, init_splits)
    rs_p = 0
    rs_n = 0
    rs_init = length(init_splits) / length(true_splits)
    for rs in resolved_splits
        if in(rs, true_splits, usemask=false)
            rs_p += 1
        else
            rs_n += 1
        end
    end
    return Dict(:rs_p => rs_p / length(true_splits), :rs_n=> rs_n / length(true_splits),
                :rs_init => rs_init)

end

"""
    _eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:yule, K=2)

Evaluate the inference of MCCs by function `f`. ARG simulation made with `(N,n0,r)`.
"""
function _eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:yule, K=2)

    arg = ARGTools.SimulateARG.simulate(N, r, n0; simtype, K)

    # Real MCCs
    rMCC = ARGTools.MCCs_pairs_from_arg(arg)[1,2];
    trees = ARGTools.trees_from_ARG(arg)

    # Inferred MCCs
    iMCC = f(arg)
    t1 = trees[1]

    out = _eval_mcc_inf(rMCC, iMCC[1,2], t1)

    out[:N] = N
    out[:n0] = n0
    out[:r] = r
    return out
end




"""
    df_fields()

Fields returned by `_eval_mcc_inf`:
- `N`, `n0`, `r`
- `ν(n)c_p(n)`: Fraction of inferred (non-)common branches, true (p) or false (n)
- `τ(n)c_p(n)`: Length of inferred (non-)common branches, true (p) or false (n) (normalized)
- `mMCCsize`: mean inferred MCC size
- `nMCC`: number of inferred MCCs
"""
function df_fields()
    return (N = Int64[], n0 = Int64[], r = Float64[], # Simulation parameters
        νc_p = Float64[], νnc_p = Float64[], # Fraction of correctly inferred common and non-common branches
        νc_n = Float64[], νnc_n = Float64[], # Fraction of incorrectly inferred common and non-common branches
        τc_p = Float64[], τnc_p = Float64[], # Total length of correctly inferred common and non-common branches
        τc_n = Float64[], τnc_n = Float64[], # Total length of incorrectly inferred common and non-common branches
        # rs_p = Float64[], rs_n = Float64[], # Number of (in)-correctly resolved splits, scaled by total number of true splits
        # rs_init = Float64[], # Initial number of splits scaled by total number of true splits
        mMCCsize = Float64[], # Mean MCC size (inferred)
        nMCC = Int64[], # Number of MCCs (inferred)
        )
end



