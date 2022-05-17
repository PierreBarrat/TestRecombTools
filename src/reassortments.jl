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
		:nrea_real_unscaled => out_inf[5],
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


function eval_reassortments_w_confidence(
	t1::Tree, t2::Tree, iMCCs, rMCCs;
	score=:branch_likelihood, neighbours=:leaves
)
	TR = 0 # true reassortments
	FR = 0 # false reassortments
	MR = 0 # missed reassortments
	#
	inf_reassortments = [lca(t1, m).label for m in iMCCs]
	real_reassortments = [lca(t1, m).label for m in rMCCs]
	conf_scores = TreeKnit.confidence_likelihood_ratio(iMCCs, t1, t2; neighbours)
	#
	out = zeros(Float64, length(inf_reassortments), 3)
	for (i, R) in enumerate(inf_reassortments)
		out[i, 1] = in(R, real_reassortments)
		out[i, 2] = 1 - out[i,1]
		if score == :branch_likelihood
			out[i, 3] = conf_scores[i]
		elseif score == :rand
			out[i, 3] = rand()
		else
			error("Unknown score type")
		end
	end
	# Removing infinity scores - they correspond to MCCs that include the root
	idx = findall(x->in(Inf, x), eachrow(out))
	out = out[findall(i->!in(i, idx), 1:size(out,1)), :]
	out
end
function read_eval_reassortments_w_confidence(folder, Nrep; score=:branch_likelihood, neighbours=:leaves)
	out_inferred = []
	out_naive = []
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			t1s = read_simulate_trees("$(folder)/$(f)/tree1.nwk", Nrep)
			t2s = read_simulate_trees("$(folder)/$(f)/tree2.nwk", Nrep)
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)
			for (rmccs, nmccs, imccs, t1, t2) in zip(rMCCs, nMCCs, iMCCs, t1s, t2s)
				out_ = eval_reassortments_w_confidence(t1, t2, imccs, rmccs; score, neighbours)
				push!(out_inferred, out_)
			end
		end
	end
	return out_inferred
end

function average_roc_curves(dat; method = :cat)
	if method == :cat
		datcat = vcat(dat...)
		if isempty(datcat) || !in(0, datcat[:,1])
			return [missing], [missing], missing
		else
			curve = roc(datcat[:,3], datcat[:,1])
			return curve.FPR, curve.TPR, AUC(curve)
		end
	elseif method == :average
		return _average_roc_curves(dat)
	end
	error("Incorrect `method` kwarg")
end

function _average_roc_curves(dat)
	dx = 0.01
	FP_grid = 0:dx:1.
	TP_ongrid = zeros(length(FP_grid))
	Z = 0
	for (j, d) in enumerate(Iterators.filter(x->size(x,1)>1, dat))
		curve = roc(d[:,3], Bool.(d[:,1]))
		if isnothing(findfirst(isnan, curve.FPR)) &&  isnothing(findfirst(isnan, curve.TPR))
		# Discard ROC with only false/true positives
			itp = LinearInterpolation(
				Interpolations.deduplicate_knots!(curve.FPR), curve.TPR
			)
			for (i, x) in enumerate(FP_grid[1:end])
				TP_ongrid[i] += itp(x)
			end
			Z += 1
		end
	end
	TP_ongrid[end] = Z

	return FP_grid, TP_ongrid/Z, sum(TP_ongrid[i]*dx/Z for (i,x) in enumerate(FP_grid[1:end-1]))
end
