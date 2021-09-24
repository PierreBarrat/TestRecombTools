function read_eval_branches(dir; filter = Dict(), Nrep = 200, verbose=true)
	df = DataFrame()
	N = length(readdir(dir))
	for (i,f) in enumerate(readdir(dir, join=true))
		verbose && print("$i / $N")
		if !isnothing(match(r"MCCs_", f))
			args = parse_outfolder(basename(f))
			if mapreduce(x->in(args[x], filter[x]), *, keys(filter), init=true)
				dat = _read_eval_branches(f, Nrep)
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

function _read_eval_branches(folder, Nrep)
	out_inf = zeros(Union{Missing, Float64}, 10)
	out_naive = zeros(Union{Missing, Float64}, 10)
	Z = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			tree = read_simulate_trees("$(folder)/$(f)/tree1.nwk", Nrep) # one is enough
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)
			for (rmccs, nmccs, imccs, t) in zip(rMCCs, nMCCs, iMCCs, tree)
				out_inf += _eval_branches(rmccs, imccs, t)
				out_naive += _eval_branches(rmccs, nmccs, t)
				Z += 1
			end
		end
	end
	out_inf = round.(out_inf / Z, digits=4)
	out_naive = round.(out_naive / Z, digits=4)
	out = Dict(
		# Real - Nmcc
		:Nmcc_real => (out_naive[10] + out_inf[10])/2,
		# Naive
		## PPV - weighed by branch length
		:τc_p_naive => out_naive[1],
		:τnc_p_naive => out_naive[2],
		:τc_n_naive => out_naive[3],
		:τnc_n_naive => out_naive[4],
		## PPV
		:νc_p_naive => out_naive[5],
		:νnc_p_naive => out_naive[6],
		:νc_n_naive => out_naive[7],
		:νnc_n_naive => out_naive[8],
		## Nmcc
		:Nmcc_naive => out_naive[9],
		# Inferred
		## PPV - weighed by branch length
		:τc_p_inf => out_inf[1],
		:τnc_p_inf => out_inf[2],
		:τc_n_inf => out_inf[3],
		:τnc_n_inf => out_inf[4],
		## PPV
		:νc_p_inf => out_inf[5],
		:νnc_p_inf => out_inf[6],
		:νc_n_inf => out_inf[7],
		:νnc_n_inf => out_inf[8],
		## Nmcc
		:Nmcc_inf => out_inf[9],
	)
	return out
end


function _eval_branches(rMCC, iMCC, t1::Tree)

    # Common branches
    τc_p = 0.; τnc_p = 0.; τc_n = 0.; τnc_n = 0.
    νc_p = 0.; νnc_p = 0.; νc_n = 0.; νnc_n = 0.
    Zn = length(nodes(t1)) - 1
    Zt = sum(skipmissing(n.tau for n in nodes(t1)))
    for n in Iterators.filter(x->!x.isroot, values(t1.lnodes))
        if RecombTools.is_branch_in_mccs(n, iMCC) # Branch predicted to be shared with the other tree
            if RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τc_p += n.tau
                νc_p += 1.
            else
                τc_n += n.tau
                νc_n += 1.
            end
        else # Branch predicted to not be shared with the other tree
            if !RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τnc_p += n.tau
                νnc_p += 1.
            else
                τnc_n += n.tau
                νnc_n += 1.
            end
        end
    end

    νc_p /= Zn; νnc_p /= Zn; νc_n /= Zn; νnc_n /= Zn
    τc_p /= Zt; τnc_p /= Zt; τc_n /= Zt; τnc_n /= Zt

    # return Dict(:τc_p => τc_p, :τnc_p => τnc_p, :τc_n => τc_n, :τnc_n => τnc_n,
    #             :νc_p => νc_p, :νnc_p => νnc_p, :νc_n => νc_n, :νnc_n => νnc_n,
    #             :Nmcc_i => length(iMCC), :Nmcc_r => length(rMCC))
    return [
    	τc_p, τnc_p, τc_n, τnc_n,
    	νc_p, νnc_p, νc_n, νnc_n,
    	length(iMCC), length(rMCC)
    ]
end


# """
#     df_fields()
# """
# function df_fields()
#     return [
#         νc_p_naive = Float64[], νnc_p_naive = Float64[], # Fraction of correctly inferred common and non-common branches
#         νc_n_naive = Float64[], νnc_n_naive = Float64[], # Fraction of incorrectly inferred common and non-common branches
#         τc_p_naive = Float64[], τnc_p_naive = Float64[], # Total length of correctly inferred common and non-common branches
#         τc_n_naive = Float64[], τnc_n_naive = Float64[], # Total length of incorrectly inferred common and non-common branches
#         νc_p_inf = Float64[], νnc_p_inf = Float64[], # ~
#         νc_n_inf = Float64[], νnc_n_inf = Float64[], # ~
#         τc_p_inf = Float64[], τnc_p_inf = Float64[], # ~
#         τc_n_inf = Float64[], τnc_n_inf = Float64[], # ~
#         Nmcc_r # Real number of MCCs
#         Nmcc_naive # Number of naive MCCs
#         Nmcc_i # Number of inferred MCCs
# 	]
# end
