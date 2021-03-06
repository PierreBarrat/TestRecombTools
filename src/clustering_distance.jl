"""
	simulate_compute_clustering_distances(dir, sim)
"""
function simulate_compute_clustering_distances(dir, sim; filter = Dict(), Nrep=-1)
	df = DataFrame()
	N = length(readdir(dir))
	for (i,f) in enumerate(readdir(dir, join=true))
		print("$i / $N - $f")
		if !isnothing(match(r"MCCs_", basename(f)))
			args = parse_outfolder(basename(f))
			if mapreduce(x->in(args[x], filter[x]), *, keys(filter), init=true)
				dat = _compute_clustering_distances(f, sim, Nrep)
				args[:naive_to_real] = dat[1]
				args[:inferred_to_real] = dat[2]
				args[:inferred_to_naive] = dat[3]
				for x in names(df)
			       df[!,x] = convert(Vector{Union{Missing, eltype(df[!,x])}}, df[!,x])
	       		end
				append!(df, DataFrame(args))
			end
		end
		print("                                            \r")
	end

	return sort!(df, [:ρ])
end

function _compute_clustering_distances(folder, sim, Nrep)
	# args = parse_outfolder(folder)
	dat = zeros(Union{Missing,Float64}, 3)
	Z = 0
	for f in readdir(folder)
		if !isnothing(match(r"rep", f)) && length(readdir(folder * "/" * f)) >= 3
			rep = parse(Int, f[4:end])
			rMCCs = read_simulate_mccs("$(folder)/$(f)/realMCCs.dat", Nrep)
			nMCCs = read_simulate_mccs("$(folder)/$(f)/naiveMCCs.dat", Nrep)
			iMCCs = read_simulate_mccs("$(folder)/$(f)/inferredMCCs.dat", Nrep)
			for i in eachindex(rMCCs)
				dat[1] += sim(rMCCs[i], nMCCs[i])
				dat[2] += sim(rMCCs[i], iMCCs[i])
				# try
				# catch
				# 	println(length(iMCCs))
				# 	println(length(rMCCs))
				# 	println(folder)
				# end
				dat[3] += sim(nMCCs[i], iMCCs[i])
				Z += 1
			end
		end
	end
	if Z == 0
		dat = [missing, missing, missing]
	end
	return dat / Z
end




#########################################################
################## SIMILARITY MEASURES ##################
#########################################################

"""
	branch_similarity(t::Tree, MCC1, MCC2)

Return fraction of branches in `t` that are predicted to be in an mcc for both `MCC1` and
  `MCC2`.
"""
function branch_similarity(t::Tree, MCC1, MCC2)
	s = 0
	for n in nodes(t)
		if TreeKnit.is_branch_in_mccs(n, MCC1) == TreeKnit.is_branch_in_mccs(n, MCC2)
			s += 1
		end
	end
	return s / length(nodes(t))
end

"""
	varinfo_similarity(MCCs...)
"""
varinfo_similarity(t::Tree, MCCs... ; scale=true) = varinfo_similarity(MCCs...; scale)
function varinfo_similarity(MCCs...; scale=true)
	leaves = sort(vcat(first(MCCs)...))
	assignments = try
		[assignment_vector(leaves, mccs) for mccs in MCCs]
	catch err
		println(MCCs[1])
		println()
		println(MCCs[2])
		println()
		println(leaves)
		error(err)
	end
	out = 0
	Z = 0
	for i in 1:length(MCCs), j in (i+1):length(MCCs)
		out += Clustering.varinfo(assignments[i], assignments[j])
		Z += 1
	end
	if scale
		out /= log(length(leaves))
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

# function get_r(ρ, n, N, simtype::Symbol)
#     if simtype == :kingman
#         return ρ * n / N
#     elseif simtype == :yule
#         return ρ / N
#     elseif simtype == :flu
#     	return ρ * n^0.2 / N
#     else
#         @error "Unrecognized `simtype`."
#     end
# end
