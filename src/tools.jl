function make_outfolder_name(N, n, ρ, cutoff, γ, nit, resolve, preresolve)
	out = "MCCs_N$(N)_n$(n)_rho$(round(ρ, sigdigits=2))_cutoff$(cutoff)"
	out *= "_gamma$(γ)_nit$(ceil(Int, nit(n)))_resolve$(resolve)_preresolve$(preresolve)"
	return out
end

function parse_outfolder(of)
	sof = String.(split(of, '_'))
	# println("pof: $of")
	args = Dict()
	for s in sof
		if s[1] == 'N'
			args[:N] = parse(Int, s[2:end])

		elseif length(s) >= 1 && s[1] == 'n' && s[2] != 'i'
			args[:n] = parse(Int, s[2:end])

		elseif length(s) >= 3 && s[1:3] == "rho"
			args[:ρ] = parse(Float64, s[4:end])

		elseif length(s) >= 6 && s[1:6] == "cutoff"
			args[:cutoff] = parse(Float64, s[7:end])

		elseif length(s) >= 5 && s[1:5] == "gamma"
			args[:γ] = parse(Float64, s[6:end])

		elseif length(s) >= 3 && s[1:3] == "nit"
			args[:nit] = parse(Int, s[4:end])

		elseif length(s) >= 7 && s[1:7] == "resolve"
			args[:resolve] = parse(Bool, s[8:end])

		elseif length(s) >= 10 && s[1:10] == "preresolve"
			args[:preresolve] = parse(Bool, s[11:end])

		end
	end

	return args
end




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

function read_simulate_trees(file::AbstractString, Nrep = -1)
	trees = []
	nwk = ""
	rep = 0
	open(file, "r") do f
		for l in eachline(f)
			if l[1] == '#'
				nwk != "" && push!(trees, TreeTools.node2tree(parse_newick(nwk)))
				nwk = ""
				if Nrep > 0 && rep >= Nrep
					break
				end
				rep += 1
			else
				nwk *= l
			end
		end
	end
	nwk != "" && push!(trees, TreeTools.node2tree(parse_newick(nwk)))

	return trees
end

function read_simulate_mccs(file::AbstractString, Nrep = -1)
	MCCs = Vector{Vector{Vector{String}}}(undef, 0)
	mcc = Vector{Vector{String}}(undef, 0)
	rep = 0
	open(file, "r") do f
		for l in eachline(f)
			if l[1] == '#'
				mcc != [] && push!(MCCs, mcc)
				mcc = []
				if Nrep > 0 && rep >= Nrep
					break
				end
				rep += 1
			else
				push!(mcc, String.(split(l, ',')))
			end
		end
	end
	mcc != [] && push!(MCCs, mcc)

	return MCCs
end

function read_simulate_results(indir::AbstractString; filter = Dict(), Nrep = 1, verbose=true)
	df = DataFrame()
	N = length(readdir(indir))
	for (i, mccfolder) in enumerate(readdir(indir, join=true))
		verbose && print("$i / $N")
		if !isnothing(match(r"MCCs_", mccfolder))
			args = parse_outfolder(basename(mccfolder))
			if mapreduce(x->in(args[x], filter[x]), *, keys(filter), init=true)
				for repfile in readdir(mccfolder)
					if !isnothing(match(r"rep",repfile)) && length(readdir(mccfolder * "/" *repfile)) >= 3
						dat = Dict()
						dat[:t1] = read_simulate_trees("$(mccfolder)/$(repfile)/tree1.nwk", Nrep)
						dat[:t2] = read_simulate_trees("$(mccfolder)/$(repfile)/tree2.nwk", Nrep)
						dat[:rMCCs] = read_simulate_mccs("$(mccfolder)/$(repfile)/realMCCs.dat", Nrep)
						dat[:nMCCs] = read_simulate_mccs("$(mccfolder)/$(repfile)/naiveMCCs.dat", Nrep)
						dat[:iMCCs] = read_simulate_mccs("$(mccfolder)/$(repfile)/inferredMCCs.dat", Nrep)
						for (k,v) in dat
							args[k] = v
						end
						append!(df, DataFrame(args))
					end
				end
			end
		end
		print("                    \r")
	end

	return sort!(df, [:ρ])
end

function _read_simulate_results(folder::AbstractString, func::Function)
	args = parse_outfolder(basename(folder))
	dat = zeros(Float64, 3)
	Z = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			# rep = parse(Int, f[4:end])
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat")
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat")
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat")
			# println(iMCCs)
			# println(func(iMCCs))
			dat[1] += mapreduce(func, +, iMCCs; init=0)
			dat[2] += mapreduce(func, +, rMCCs; init=0)
			dat[3] += mapreduce(func, +, nMCCs; init=0)
			Z += length(rMCCs)
		end
	end
	if Z == 0
		dat = [missing, missing, missing]
	end
	return dat / Z
end

"""
	read_simulate_results(dir, func::Dict)
"""
function read_simulate_results(dir, func::Dict)
	df = DataFrame()
	N = length(readdir(dir))
	for (i,f) in enumerate(readdir(dir, join=true))
		print("$i / $N")
		if !isnothing(match(r"MCCs_", f))
			args = parse_outfolder(basename(f))
			for (name, fun) in func
				dat = _read_simulate_results(f, fun)
				args[Symbol(name, :_inferred)] = dat[1]
				args[Symbol(name, :_real)] = dat[2]
				args[Symbol(name, :_naive)] = dat[3]
			end

			for x in names(df)
		       df[!,x] = convert(Vector{Union{Missing, eltype(df[!,x])}}, df[!,x])
       		end
			append!(df, DataFrame(args))
		end
		print("                    \r")
	end

	return sort!(df, [:ρ])
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
"""
	remove_branches!(t::Tree, c::Real, N::Int)
"""
function remove_branches!(t::Tree, c::Real, N::Int)
	if c > 0
		remove_branches!(t, Exponential(c*N))
	end
	return nothing
end



function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / N
    elseif simtype == :yule
        return ρ / N
    elseif simtype == :flu
    	return ρ * n^0.2 / N
    else
        @error "Unrecognized `simtype`."
    end
end
