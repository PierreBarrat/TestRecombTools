function eval_reassortments(t::Tree, iMCCs, rMCCs)
	TR = 0 # true reassortments
	FR = 0 # false reassortments
	MR = 0 # missed reassortments
	#
	inf_reassortments = [lca(t, m).label for m in iMCCs]
	real_reassortments = [lca(t, m).label for m in rMCCs]
	# Removing root: uninteresting case
	root = t.root.label
	deleteat!(inf_reassortments, findall(==(root), inf_reassortments))
	deleteat!(real_reassortments, findall(==(root), real_reassortments))
	#
	TR = length(intersect(inf_reassortments, real_reassortments))
	FR = length(setdiff(inf_reassortments, real_reassortments))
	MR = length(setdiff(real_reassortments, inf_reassortments))
	#
	return [TR, FR, MR, length(inf_reassortments), length(real_reassortments)]
end

function eval_reassortment(t::Tree, iMCC::Vector{<:AbstractString}, rMCCs)
	inf_reassortment = lca(t, iMCC).label
	real_reassortments = [lca(t, m).label for m in rMCCs]
	return in(inf_reassortment, real_reassortments)
end

function _read_eval_reassortments(folder, Nrep)
	out_inf = zeros(Union{Missing, Float64}, 5)
	out_naive = zeros(Union{Missing, Float64}, 5)
	Z_FR_naive = 0; Z_FR_inf = 0
	Z_MR_naive = 0; Z_MR_inf = 0
	Z = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			tree = read_simulate_trees("$(folder)/$(f)/tree1.nwk", Nrep) # one is enough
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)
			for (rmccs, nmccs, imccs, t) in zip(rMCCs, nMCCs, iMCCs, tree)
				out = eval_reassortments(t, imccs, rmccs)
				out_inf += out
				out[4] > 0 && (Z_FR_inf += 1) # FR --> only if we inferred a rea
				out[5] > 0 && (Z_MR_inf += 1) # MR --> only if there is a rea to miss
				#
				out = eval_reassortments(t, nmccs, rmccs)
				out_naive += out
				out[4] > 0 && (Z_FR_naive += 1) # FR --> only if we inferred a rea
				out[5] > 0 && (Z_MR_naive += 1) # MR --> only if there is a rea to miss

				Z += 1
			end
		end
	end

	out = Dict(
		:TR_naive => out_naive[1], #/ Z_FR_naive,
		:FR_naive => out_naive[2], #/ Z_FR_naive,
		:MR_naive => out_naive[3], #/ Z_MR_naive,
		:nrea_naive => out_naive[4] / Z,
		#
		:TR_inf => out_inf[1], #/ Z_FR_inf,
		:FR_inf => out_inf[2], #/ Z_FR_inf,
		:MR_inf => out_inf[3], #/ Z_MR_inf,
		:nrea_inf => out_inf[4] / Z,
		#
		:nrea_real => out_inf[5] / Z,
	)
	return out
end

function read_eval_reassortments(dir; filter = Dict(), Nrep = 200, verbose=true)
	df = DataFrame()
	N = length(readdir(dir))
	for (i,f) in enumerate(readdir(dir, join=true))
		verbose && print("$i / $N")
		if !isnothing(match(r"MCCs_", f))
			args = parse_outfolder(basename(f))
			if mapreduce(x->in(args[x], filter[x]), *, keys(filter), init=true)
				dat = _read_eval_reassortments(f, Nrep)
				for (k,v) in dat
					args[k] = v
				end
				append!(df, DataFrame(args))
			end
		end
		print("                    \r")
	end

	return sort!(df, [:ρ])
end
