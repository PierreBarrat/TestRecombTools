function _read_simulate_splits(folder::AbstractString, Nrep=200)
	args = parse_outfolder(basename(folder))
	TP = zeros(Float64, 3)
	power = zeros(Float64, 3)
	Z = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			trees1 = read_simulate_trees("$(folder)/$(f)/tree1.nwk", Nrep)
			trees2 = read_simulate_trees("$(folder)/$(f)/tree2.nwk", Nrep)
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)

			#
			for (rep, (t1, t2, iM, rM, nM)) in enumerate(zip(trees1, trees2, iMCCs, rMCCs, nMCCs))
				nS_i = RecombTools.new_splits(t1, iM, t2)
				nS_r = RecombTools.new_splits(t1, rM, t2)
				nS_n = RecombTools.new_splits(t1, nM, t2)

				TP[1] += count(in(nS_r), nS_i) / length(nS_i)
				TP[2] += count(in(nS_r), nS_r) / length(nS_r)
				TP[3] += count(in(nS_r), nS_n) / length(nS_n)

				missing_splits = 2*length(leaves(t1)) - 1 - length(nodes(t1))
				power[1] += count(in(nS_r), nS_i) / missing_splits
				power[2] += count(in(nS_r), nS_r) / missing_splits
				power[3] += count(in(nS_r), nS_n) / missing_splits
				Z += 1
				rep > Nrep && break
			end
			# Z += length(rMCCs)
		end
	end
	if Z == 0
		TP = [missing, missing, missing]
	end
	return TP / Z, power / Z
end

"""

"""
function read_simulate_splits(dir; filter = Dict(), Nrep = 200, verbose=true)
	df = DataFrame()
	N = length(readdir(dir))
	for (i,f) in enumerate(readdir(dir, join=true))
		verbose && print("$i / $N")
		if !isnothing(match(r"MCCs_", f))
			args = parse_outfolder(basename(f))
			if mapreduce(x->in(args[x], filter[x]), *, keys(filter), init=true)
				dat = _read_simulate_splits(f, Nrep)
				# TP
				args[:TP_splits_inferred] = dat[1][1]
				args[:TP_splits_real] = dat[1][2]
				args[:TP_splits_naive] = dat[1][3]
				# power
				args[:power_splits_inferred] = dat[2][1]
				args[:power_splits_real] = dat[2][2]
				args[:power_splits_naive] = dat[2][3]
				for x in names(df)
			       df[!,x] = convert(Vector{Union{Missing, eltype(df[!,x])}}, df[!,x])
	       		end
				append!(df, DataFrame(args))
			end
		end
		print("                    \r")
	end

	return sort!(df, [:ρ])
end
