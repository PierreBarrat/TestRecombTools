function _read_simulate_splits(folder::AbstractString, Nrep=200)
	# args = parse_outfolder(basename(folder))
	TP = zeros(Union{Missing,Float64}, 3)
	power = zeros(Union{Missing, Float64}, 3)
	power_neg = zeros(Union{Missing, Float64}, 2)
	Z_ppv = zeros(Int, 3)
	Z_power = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			trees1 = read_simulate_trees("$(folder)/$(f)/tree1.nwk", Nrep)
			trees2 = read_simulate_trees("$(folder)/$(f)/tree2.nwk", Nrep)
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)

			#
			for (rep, (t1, t2, iM, rM, nM)) in enumerate(zip(trees1, trees2, iMCCs, rMCCs, nMCCs))
				nS_i = TreeKnit.new_splits(t1, iM, t2)
				nS_r = TreeKnit.new_splits(t1, rM, t2)
				nS_n = TreeKnit.new_splits(t1, nM, t2)

				TP_i = count(in(nS_r), nS_i)
				TP_r = count(in(nS_r), nS_r)
				TP_n = count(in(nS_r), nS_n)

				FP_i = length(nS_i) - TP_i
				FP_n = length(nS_n) - TP_n

				if length(nS_i) != 0
					TP[1] += TP_i / length(nS_i)
					Z_ppv[1] += 1
				end
				if length(nS_r) != 0
					TP[2] += TP_r / length(nS_r)
					Z_ppv[2] += 1
				end
				if length(nS_n) != 0
					TP[3] += TP_n / length(nS_n)
					Z_ppv[3] += 1
				end

				missing_splits = 2*length(leaves(t1)) - 1 - length(nodes(t1))
				power[1] += TP_i / missing_splits
				power[2] += TP_r / missing_splits
				power[3] += TP_n / missing_splits

				power_neg[1] += FP_i / missing_splits
				power_neg[2] += FP_n / missing_splits
				Z_power += 1
				rep > Nrep && break
			end
		end
	end
	for (i,z) in enumerate(Z_ppv)
		if z < Nrep / 5
			TP[i] = missing
			power[i] = missing
		end
	end
	return TP ./ Z_ppv, power / Z_power, power_neg / Z_power
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
				#
				args[:power_neg_splits_inferred] = dat[3][1]
				args[:power_neg_splits_naive] = dat[3][2]
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
